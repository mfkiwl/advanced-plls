function [received_signal, los_phase, psi, ps_realization] = ...
    get_received_signal(C_over_N0_dBHz,S4,tau0,simulation_time, ...
    settling_time,scint_model, is_refractive_effects_removed)
% get_received_signal
% Simulates the baseband received signal, including ionospheric 
% scintillation effects, thermal noise, and line-of-sight (LOS) phase 
% dynamics, for a given scintillation model.
%
% Syntax:
%   [received_signal] = get_received_signal(C_over_N0_dBHz, S4,
%       tau0, simulation_time, settling_time, scint_model, is_refractive_effects_removed)
%   [received_signal, phi_LOS] = get_received_signal(C_over_N0_dBHz, S4,
%       tau0, simulation_time, settling_time, scint_model, is_refractive_effects_removed)
%   [received_signal, phi_LOS, psi] = ...
%       get_received_signal(C_over_N0_dBHz, S4,tau0, simulation_time, ...
%       settling_time, scint_model)
%   [received_signal, phi_LOS, psi, ps_realization] = ...
%       get_received_signal(C_over_N0_dBHz, S4,tau0, simulation_time, ...
%       settling_time, scint_model, is_refractive_effects_removed)
%
% Description:
%   This function generates a baseband received signal by combining:
%   - Line-of-sight (LOS) phase dynamics, modeled with Doppler shift 
%     and drift.
%   - Thermal noise based on the carrier-to-noise density ratio (C/N₀).
%   - Ionospheric scintillation effects modeled using either the Cornell 
%     Scintillation Model (CSM) or the Multi-Frequency Phase Screen Model 
%     (MFPSM).
%
%   The additional outputs, `psi` and `ps_realization`, are dependent on 
%   the selected scintillation model:
%   - When `CSM` is chosen, only the scintillation field (`psi`) is 
%     returned.
%   - When `MFPSM` is chosen, both the scintillation field (`psi`) and the 
%     phase screen realization (`ps_realization`) are returned.
%   - When 'none' is chosen, the signal is generated without ionospheric 
%     scintillation effects. Only thermal noise and LOS phase dynamics 
%     are applied.
%
% Inputs:
%   C_over_N0_dBHz   - Carrier-to-noise density ratio (C/N₀) in dB-Hz.
%   S4               - Scintillation index (0 <= S4 <= 1), representing 
%                      the severity of amplitude fluctuations caused by 
%                      ionospheric scintillation.
%   tau0             - Signal intensity decorrelation time in seconds, 
%                      describing temporal scintillation variability.
%   simulation_time  - Total simulation time in seconds.
%   settling_time    - Amount of time in seconds during which the receiver 
%                      is not subjected to ionospheric scintillation. 
%                      Used to ensure convergence of the Autoregressive 
%                      Kalman Filter before introducing scintillation.
%   scint_model      - Scintillation model to use. Must be either:
%                      - 'CSM': Cornell Scintillation Model.
%                      - 'MFPSM': Multi-Frequency Phase Screen Model.
%                      - 'none': No scintillation effects are applied.
%   is_refractive_effects_removed (Optional)
%                    - Logical flag to toggle the removal of 
%                      refractive effects (true/false). Defaults to `false` 
%                      if not provided. Not applicable for 
%                      scint_model = 'CSM' or 'none'.
% Outputs:
%   received_signal  - Baseband received signal (complex-valued), 
%                      incorporating ionospheric scintillation effects, 
%                      thermal noise, and LOS phase dynamics.
%   los_phase          - Line-of-sight phase time series (radians) for the 
%                      simulated signal.
%   psi              - Complex scintillation field. For 'CSM', it contains 
%                      only diffractive effects. For 'MFPSM', it may also 
%                      include refractive effects depending on the value of 
%                      `is_refractive_effects_removed`.
%   ps_realization   - Phase screen realization. Only meaningful for 
%                      'MFPSM'; empty for 'CSM'.
%
% Notes:
%   - Some parameters are fixed:
%       - Integration time (T_I): 0.01 seconds.
%       - Doppler frequency shift (fd): 1000 Hz.
%       - Doppler drift rate (fdr): 0.94 Hz/s.
%   - The receiver mean signal power is set to 1.
%   - Thermal noise is computed based on receiver bandwidth (B = 20 MHz).
%   - Currently, this function supports single-frequency simulation only. 
%     Future updates will extend it to handle multi-frequency scenarios.
%   - The settling period (`settling_time`) is used to ensure convergence of the
%     signal dynamics before introducing scintillation effects. For `'none'`, 
%     this period will still produce thermal noise and LOS phase dynamics without 
%     scintillation effects.
%
% Example 1:
%   % Simulate a received signal with:
%   % C/N₀ = 45 dB-Hz, S4 = 0.8, tau0 = 0.7 seconds, simulation time of 
%   % 300 seconds, settling time of 50 seconds, using the CSM scintillation 
%   % model:
%   [received_signal, phi_LOS, psi] = get_received_signal(45, 0.8, 0.7,
%       300, 50, 'CSM');
%
% Example 2:
%   % Simulate a received signal with:
%   % C/N₀ = 45 dB-Hz, S4 = 0.8, tau0 = 0.7 seconds, simulation time of 
%   % 300 seconds, settling time of 50 seconds, using the MFPSM 
%   % scintillation model:
%   [received_signal, phi_LOS, psi, ps_realization] = ...
%       get_received_signal(45, 0.8, 0.7, 300, 50, 'MFPSM');
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

% Input Validation
validateattributes(C_over_N0_dBHz, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite', 'nonnan', 'nonempty'}, ...
    'get_received_signal', 'C_over_N0_dBHz');
validateattributes(S4, {'numeric'}, {'scalar', '>=', 0, '<=', 1}, ...
    'get_received_signal', 'S4');
validateattributes(tau0, {'numeric'}, {'scalar', 'real', 'positive', ...
    'finite', 'nonnan'}, 'get_received_signal', 'tau0');
validateattributes(simulation_time, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite', 'nonnan', 'nonempty'}, ...
    'get_received_signal', 'simulation_time');
validateattributes(settling_time, {'numeric'}, {'scalar', 'real', ...
    'positive', 'finite', 'nonnan'}, 'get_received_signal', 'settling_time');

% Ensure settling time does not exceed simulation time
if settling_time > simulation_time
    error('get_received_signal:InvalidInput', ...
        'Settling time (%g) must not exceed simulation time (%g).', ...
        settling_time, simulation_time);
end

% Validate scintillation model
validateattributes(scint_model, {'char', 'string'}, {'nonempty'}, ...
    'get_received_signal', 'scint_model');
if ~ismember(scint_model, {'CSM', 'MFPSM', 'none'})
    error('get_received_signal:InvalidScintModel', ['Invalid scintillation ' ...
        'model. Must be ''CSM'', ''MFPSM'', or ''none''.']);
end

% Validate is_refractive_effects_removed (optional, default = false)
if nargin >= 7
    validateattributes(is_refractive_effects_removed, {'logical'}, {'scalar','nonnan'}, ...
    'get_received_signal', 'is_refractive_effects_removed');
else
    is_refractive_effects_removed = false;
end

if ismember(scint_model, {'CSM', 'none'}) && is_refractive_effects_removed == true
    warning('get_received_signal:UnusedFlag', ...
        ['The flag "is_refractive_effects_removed" is not applicable ' ...
        'for the selected scintillation model ("%s").'], scint_model);
end

% Fixed Parameters
% Sampling time of the prompt correlator signal after the integrate and
% dump.
T_I = 0.01;
% Initial phase of the line-of-sight (LOS) signal dynamics.
los_phase_0 = 0;
% Doppler frequency shift in Hz.
fd = 0;
% Doppler frequency drift rate in Hz/s
fdr = 0;
% Mean power of the received signal
rx_mean_power = 1;
% Receiver Bandwidth
B = 2e7;

% Get the values of the line-of-sight phase shift for a given simulation 
% time, sampling time, initial phase, Doppler frequency shift and Doppler 
% frequency drift rate.
los_phase = get_los_phase(simulation_time,T_I,los_phase_0,fd,fdr);
% Generate discrete thermal noise based on simulation time, sampling time, 
% receiver mean power, carrier-to-noise ratio, and receiver bandwidth.
thermal_noise = get_thermal_noise(simulation_time,T_I,rx_mean_power, ...
    C_over_N0_dBHz,B);

% Select the scintillation model for generating the received baseband 
% signal:
% - 'CSM': Cornell Scintillation Model for ionospheric scintillation 
%   effects.
% - 'MFPSM': Multi-Frequency Phase Screen Model for ionospheric 
%   scintillation effects.
% If an invalid model is specified, the function raises an error and 
% exits.
% TODO: Currently, this function supports single-frequency simulation 
% only. Future updates will extend it to handle multi-frequency 
% scenarios.
switch scint_model
    case 'CSM'
        % Generate the ionospheric scintillation effects using the 
        % Cornell Scintillation Model (CSM). This model captures only 
        % the diffractive effects of ionospheric scintillation.
        psi = get_csm_data(S4,tau0,simulation_time,T_I);
        % The refractive effects, modeled in the MFPSM as phase screen 
        % realizations, are not included in the CSM and are set to zero.
        ps_realization = [];
    case 'MFPSM'
        % Get a complex field that represents the ionospheric scintillation
        % effects using the multi-frequency phase screen model, as well as
        % the phase screen realization used to generate the aforementioned
        % complex field.
        [psi, ps_realization] = get_mfpsm_data(S4,tau0, ...
            simulation_time,T_I);
        % Remove refractive effects if `is_refractive_effects_removed` 
        % is true by applying the negative phase screen realization
        % to `psi`.
        if is_refractive_effects_removed
            psi = psi .* exp(-1j * ps_realization);
        end
     case 'none'
        % No scintillation effects
        psi = ones(simulation_time / T_I, 1);
        ps_realization = [];
end

% Add an initial settling period with no ionospheric scintillation effects.
psi_settled = zeros(simulation_time / T_I, 1);
psi_settled(1:settling_time / T_I) = 1 + 0j;
psi_settled(settling_time / T_I + 1:end) = ...
    psi(settling_time / T_I + 1:end);

% Construct the baseband received signal
received_signal = sqrt(rx_mean_power)*psi_settled .* exp(1j*los_phase) + ...
    thermal_noise;
end
function [received_signal, los_phase, psi_settled, diffractive_phase_settled, refractive_phase_settled] = get_received_signal(C_over_N0_dBHz, scint_model, doppler_profile, varargin)
% get_received_signal
% Simulates the baseband received signal, including ionospheric scintillation
% effects, thermal noise, and line-of-sight (LOS) phase dynamics.
%
% This function generates a baseband received signal by combining:
%   - LOS phase dynamics (via Doppler shift and drift),
%   - Thermal noise computed from the carrier-to-noise density ratio (C/N₀),
%   - Ionospheric scintillation effects modeled using either:
%       * the Cornell Scintillation Model (CSM) [requires 'S4' and 'tau0'],
%       * the Two-Component Power Law Model (TPPSM) [requires 'tppsm_scenario'],
%       * or no scintillation effects ('none').
%
% Optional Name-Value Pair Inputs:
%   'simulation_time'              - Duration of simulation in seconds (default: 300 s)
%   'settling_time'                - Settling time in seconds (default: 50 s)
%   'sampling_interval'            - Integration time in seconds (default: 0.01 s)
%   'is_refractive_effects_removed'- Logical flag (default: true)
%   'S4'                           - Scintillation index (0<=S4<=1); required if scint_model is 'CSM'
%   'tau0'                         - Signal intensity decorrelation time (s); required if scint_model is 'CSM'
%   'tppsm_scenario'               - TPPSM scenario ('weak', 'moderate', or 'strong'); required if scint_model is 'TPPSM'
%   'is_enable_cmd_print'          - logical value that configures whether command
%                                    lines would appear on the regarding of the 
%                                    usage of external rhof_veff_ratio.
%
% Inputs:
%   C_over_N0_dBHz  - Carrier-to-noise density ratio in dB-Hz.
%   scint_model     - Scintillation model to use ('CSM', 'TPPSM', or 'none').
%   doppler_profile - Array of coefficients for the Taylor series expansion 
%                     of LOS phase dynamics.
%
% Outputs:
%   received_signal - Baseband received signal (complex-valued).
%   los_phase       - LOS phase time series (radians).
%   psi             - Complex scintillation field. For 'CSM', contains only diffractive
%                     effects; for 'TPPSM', contains the model's combined effects.
%   refractive_phase  - Phase screen realization (only meaningful for TPPSM; empty for CSM or 'none').
%   diffractive_phase - Diffractive phase of the ionospheric scintillation
%                       complex time series
%
% Example 1 (using CSM):
%   [rx_sig, los_phase, psi] = get_received_signal(45, 'CSM', [0,1000,0.94,0.01], ...
%         'S4', 0.8, 'tau0', 0.7, 'simulation_time', 300, 'settling_time', 50, 'is_refractive_effects_removed', false);
%
% Example 2 (using TPPSM):
%   [rx_sig, los_phase, psi, refractive_phase] = get_received_signal(45, 'TPPSM', [0,1000,0.94,0.01], ...
%         'tppsm_scenario', 'Moderate', 'simulation_time', 300, 'settling_time', 50, 'is_refractive_effects_removed', true);
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

% Input Parsing
p = inputParser;
addRequired(p, 'C_over_N0_dBHz', @(x) isnumeric(x) && isscalar(x) && (x > 0));
addRequired(p, 'scint_model', @(x) ischar(x) || isstring(x));
addRequired(p, 'doppler_profile', @(x) isnumeric(x) && ~isempty(x));

% Optional name-value parameters
addParameter(p, 'simulation_time', 300, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p, 'settling_time', 50, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p, 'sampling_interval', 0.01, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p, 'is_refractive_effects_removed', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'S4', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 1)));
addParameter(p, 'tau0', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && (x > 0)));
addParameter(p, 'tppsm_scenario', [], @(x) isempty(x) || ismember(x, {'weak', 'moderate', 'strong'}));
addParameter(p, 'is_enable_cmd_print', true, @(x) validateattributes(x, {'logical'}, {'nonempty'}));

parse(p, C_over_N0_dBHz, scint_model, doppler_profile, varargin{:});

simulation_time = p.Results.simulation_time;
settling_time = p.Results.settling_time;
sampling_interval = p.Results.sampling_interval;
is_refractive_effects_removed = p.Results.is_refractive_effects_removed;
S4 = p.Results.S4;
tau0 = p.Results.tau0;
tppsm_scenario = p.Results.tppsm_scenario;
is_enable_cmd_print = p.Results.is_enable_cmd_print;

scint_model = validatestring(p.Results.scint_model, {'CSM', 'TPPSM', 'none'}, mfilename, 'scint_model');

% Validate required parameters based on selected scintillation model
switch scint_model
    case 'CSM'
        if isempty(S4) || isempty(tau0)
            error('get_received_signal:MissingParameters', ...
                'For scint_model ''CSM'', both ''S4'' and ''tau0'' must be provided.');
        end
    case 'TPPSM'
        if isempty(tppsm_scenario)
            error('get_received_signal:MissingParameters', ...
                'For scint_model ''TPPSM'', the ''tppsm_scenario'' parameter must be provided.');
        end
    case 'none'
        % No additional parameters required.
end

% Validate that settling_time does not exceed simulation_time
if settling_time > simulation_time
    error('get_received_signal:InvalidInput', ...
        'Settling time (%g s) must not exceed simulation time (%g s).', settling_time, simulation_time);
end

% Fixed parameters
rx_mean_power = 1;
sampling_frequency_before_correlation = 2e7;
B = sampling_frequency_before_correlation/2; % Receiver bandwidth in Hz

% Generate LOS phase using the provided doppler_profile
los_phase = get_los_phase(simulation_time, sampling_interval, doppler_profile);

% Generate thermal noise
thermal_noise = get_thermal_noise(simulation_time, sampling_interval, rx_mean_power, C_over_N0_dBHz, B);

% Select scintillation model for baseband signal generation
switch scint_model
    case 'TPPSM'
        % Call get_tppsm_data with the provided tppsm_scenario.
        [psi, refractive_phase] = get_tppsm_data(tppsm_scenario, 'is_enable_cmd_print', is_enable_cmd_print, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 0.5);
        diffractive_phase = wrapToPi(unwrap(angle(psi)) - refractive_phase);
        % Remove refractive effects if requested.
        if is_refractive_effects_removed
            psi = psi .* exp(-1j * refractive_phase);
        end
        
    case 'CSM'
        % Use the Cornell Scintillation Model with S4 and tau0.
        csm_params = struct('S4', S4, 'tau0', tau0, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval);
        psi = get_csm_data(csm_params);
        refractive_phase = [];
        diffractive_phase = angle(psi);
    case 'none'
        psi = ones(round(simulation_time / sampling_interval), 1);
        refractive_phase = [];
        diffractive_phase = zeros(size(psi));
end

% Apply settling period: initial period without scintillation effects
nSamples = round(simulation_time / sampling_interval);
psi_settled = zeros(nSamples, 1);
settling_samples = round(settling_time / sampling_interval);

diffractive_phase_settled = zeros(nSamples, 1);
refractive_phase_settled = zeros(nSamples, 1);

psi_settled(1:settling_samples) = 1 + 0j;
psi_settled(settling_samples+1:end) = psi(settling_samples+1:end);

diffractive_phase_settled(1:settling_samples) = 0;
diffractive_phase_settled(settling_samples+1:end) = diffractive_phase(settling_samples+1:end);

if ~isempty(refractive_phase)
    refractive_phase_settled(1:settling_samples) = 0;
    refractive_phase_settled(settling_samples+1:end) = refractive_phase(settling_samples+1:end);
end

% Construct the baseband received signal
received_signal = sqrt(rx_mean_power) * psi_settled .* exp(1j * los_phase) + thermal_noise;
end

% Script: script_KFAR_RMSE_monte_carlo_with_process_noise_sweep.m
%
% Description:
%   This script performs a Monte Carlo simulation where the process noise variance is
%   swept over a logarithmically spaced range. For each value of process noise variance and 
%   for a given number of Monte Carlo runs, the script:
%     1. Generates received signal data for both CSM and TPPSM scenarios under severe scintillation.
%     2. Configures the Kalman PLL (KFAR) and standard KF algorithms with different adaptive 
%        update options.
%     3. Computes the corresponding state estimates.
%     4. Computes the RMSE metrics for three phase types:
%          - LOS phase error,
%          - Scintillation phase error (only for AR-augmented estimates),
%          - Joint phase error (LOS + scintillation).
%     5. Stores these RMSE values for each algorithm variant:
%          - KF-AR      : KF with AR augmentation (no adaptive update)
%          - AKF-AR     : KF with AR augmentation using simplified adaptation (hard_limited = false)
%          - AHL-KF-AR  : KF with AR augmentation using simplified adaptation (hard_limited = true)
%          - KF-std     : Standard KF (no AR augmentation; training_scint_model = 'none')
%          - AKF-std    : Standard KF using simplified adaptation (hard_limited = false)
%          - AHL-KF-std : Standard KF using simplified adaptation (hard_limited = true)
%
%   The RMSE is computed over the valid time vector (i.e. after a settling time) by comparing
%   the estimated phase with the “true” phase (for LOS: unwrap(angle(rx_sig)) - los_phase; for scintillation: angle(psi); and for joint: their sum).
%
%   All RMSE values are stored in a structured array named "results" that organizes the metrics by
%   scenario (CSM, TPPSM), phase type (los, scint, joint), and algorithm variant.
%
% [1] R. A. M. Lopes, F. Antreich, F. Fohlmeister, M. Kriegel and H. K. Kuga, "Ionospheric 
%     Scintillation Mitigation With Kalman PLLs Employing Radial Basis Function Networks," 
%     in IEEE Transactions on Aerospace and Electronic Systems, vol. 59, no. 5, pp. 6878-6893,
%     Oct. 2023, doi: 10.1109/TAES.2023.3281431
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', 'libs')));

%% Overall Parameters
% Simulation Parameters
doppler_profile      = [0, 1000, 0.94];  % Doppler profile parameters
L1_C_over_N0_dBHz    = 42;               % Carrier-to-noise ratio (in dBHz)
simulation_time      = 300;              % Total simulation time (sec)
S4                   = 0.8;              % S4 index for CSM
tau0                 = 0.5;              % tau_0 for CSM
settling_time        = 50;               % Settling time (sec) -> Amount 
                                         % of time without scintillation 
                                         % happening in the beginning of the simulation.

% Sampling Parameters
sampling_interval    = 0.01;             % Sampling interval (sec) – 100 Hz sampling rate
time_vector          = (sampling_interval:sampling_interval:simulation_time);
valid_time_vector    = ((settling_time/sampling_interval + 1) : simulation_time/sampling_interval).';

% AR Model Parameters
var_model_order        = 5;              % AR model order
initial_process_noise_variance = 2.6e-6; % Process noise variance (e.g., from reference [1])


% Initial state distribution boundaries. It is assumed that the initial
% values are sampled from a uniform distribution.
initial_states_distributions_boundaries = { {[-pi, pi], [-5, 5], [-0.1, 0.1]} };

% Training Data Configuration
training_simulation_time = 1500;         % Training simulation time (sec)
training_data_config_csm = struct(...
    'scintillation_model', 'CSM', ...
    'S4', S4, ...
    'tau0', tau0, ...
    'simulation_time', training_simulation_time, ...
    'sampling_interval', sampling_interval);
training_data_config_tppsm = struct(...
    'scintillation_model', 'TPPSM', ...
    'scenario', 'Severe', ...
    'simulation_time', training_simulation_time, ...
    'is_refractive_effects_removed', true, ...
    'sampling_interval', sampling_interval);
training_data_config_none = struct(...
    'scintillation_model', 'none', ...
    'sampling_interval', sampling_interval);

% Miscellaneous Settings
cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');
scint_index = length(doppler_profile) + 1;  % Column index for scintillation phase (LOS is in column 1)
is_enable_cmd_print_get_kf_cfg_fnc = false;

%% Define adaptive configuration structures
% For the simplified adaptive update, we require L1_C_over_N0_dBHz, sampling_interval, and threshold.
hard_limited_threshold = 38;
% For AR (KFAR) estimates:
adaptive_config_KF_AR = struct('algorithm', 'none', 'hard_limited', false);
adaptive_config_AKF_AR = struct('algorithm', 'simplified', 'hard_limited', false, ...
                                'L1_C_over_N0_dBHz', L1_C_over_N0_dBHz, ...
                                'sampling_interval', sampling_interval, ...
                                'threshold', hard_limited_threshold);
adaptive_config_AHL_KF_AR = struct('algorithm', 'simplified', 'hard_limited', true, ...
                                   'L1_C_over_N0_dBHz', L1_C_over_N0_dBHz, ...
                                   'sampling_interval', sampling_interval, ...
                                   'threshold', hard_limited_threshold);
% For standard KF estimates (training_scint_model = 'none'):
adaptive_config_KF_std = struct('algorithm', 'none', 'hard_limited', false);
adaptive_config_AKF_std = struct('algorithm', 'simplified', 'hard_limited', false, ...
                                 'L1_C_over_N0_dBHz', L1_C_over_N0_dBHz, ...
                                 'sampling_interval', sampling_interval, ...
                                 'threshold', hard_limited_threshold);
adaptive_config_AHL_KF_std = struct('algorithm', 'simplified', 'hard_limited', true, ...
                                    'L1_C_over_N0_dBHz', L1_C_over_N0_dBHz, ...
                                    'sampling_interval', sampling_interval, ...
                                    'threshold', hard_limited_threshold);

%% Initializing the general_config structs for `CSM`, `TPPSM` and `none`.

discrete_wiener_model_config = {1, 3, 0.01, [0, 0, 1e-2], 1};
general_config_csm = struct( ...
  'discrete_wiener_model_config', {discrete_wiener_model_config}, ...
  'scintillation_training_data_config', training_data_config_csm, ...
  'var_minimum_order', var_model_order, ...
  'var_maximum_order', var_model_order, ...
  'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
  'initial_states_distributions_boundaries', { {[-pi, pi], [-5, 5], [-0.1, 0.1]} }, ...
  'real_doppler_profile', doppler_profile, ...
  'is_use_cached_settings', false, ...
  'is_generate_random_initial_estimates', true ...
);

general_config_tppsm = general_config_csm;
general_config_tppsm.scintillation_training_data_config = training_data_config_tppsm;

general_config_none = general_config_csm;
general_config_none.scintillation_training_data_config = training_data_config_none;

%% Main Monte Carlo Loop

monte_carlo_runs = 2;
pnv_amount = 2;
pnv_array = logspace(-4,-9,pnv_amount);
% Preallocate a results structure for RMSE metrics
results = struct();

% For CSM scenario:
results.csm.los.kf_ar           = zeros(monte_carlo_runs,pnv_amount);
results.csm.los.akf_ar          = zeros(monte_carlo_runs,pnv_amount);
results.csm.los.ahl_kf_ar       = zeros(monte_carlo_runs,pnv_amount);
results.csm.scint.kf_ar         = zeros(monte_carlo_runs,pnv_amount);
results.csm.scint.akf_ar        = zeros(monte_carlo_runs,pnv_amount);
results.csm.scint.ahl_kf_ar     = zeros(monte_carlo_runs,pnv_amount);
results.csm.joint.kf_ar         = zeros(monte_carlo_runs,pnv_amount);
results.csm.joint.akf_ar        = zeros(monte_carlo_runs,pnv_amount);
results.csm.joint.ahl_kf_ar     = zeros(monte_carlo_runs,pnv_amount);
results.csm.joint.kf_std        = zeros(monte_carlo_runs,pnv_amount);
results.csm.joint.akf_std       = zeros(monte_carlo_runs,pnv_amount);
results.csm.joint.ahl_kf_std    = zeros(monte_carlo_runs,pnv_amount);

% For TPPSM scenario:
results.tppsm.los.kf_ar         = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.los.akf_ar        = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.los.ahl_kf_ar     = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.scint.kf_ar       = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.scint.akf_ar      = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.scint.ahl_kf_ar   = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.joint.kf_ar       = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.joint.akf_ar      = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.joint.ahl_kf_ar   = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.joint.kf_std      = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.joint.akf_std     = zeros(monte_carlo_runs,pnv_amount);
results.tppsm.joint.ahl_kf_std  = zeros(monte_carlo_runs,pnv_amount);

%% Monte Carlo Loop 
for pnv_idx = 1:pnv_amount
    fprintf('Processing noise variance %d/%d: %g\n', pnv_idx, pnv_amount, pnv_array(pnv_idx));

    discrete_wiener_model_config{1,4}(3) = pnv_array(pnv_idx);
    
    general_config_csm.discrete_wiener_model_config = discrete_wiener_model_config;
    general_config_tppsm.discrete_wiener_model_config = discrete_wiener_model_config;
    general_config_none.discrete_wiener_model_config = discrete_wiener_model_config;
    
    [~, ~] = get_kalman_pll_config(general_config_csm, cache_dir, is_enable_cmd_print_get_kf_cfg_fnc);
    [~, ~] = get_kalman_pll_config(general_config_tppsm, cache_dir, is_enable_cmd_print_get_kf_cfg_fnc);
    [kf_cfg, ~] = get_kalman_pll_config(general_config_none, cache_dir, is_enable_cmd_print_get_kf_cfg_fnc);
    
    for m = 1:monte_carlo_runs
        fprintf('  Monte Carlo run: %d/%d\n', m, monte_carlo_runs);

        rng(m);

        init_estimates_csm = get_initial_estimates(general_config_csm, kf_cfg);
        init_estimates_tppsm = get_initial_estimates(general_config_tppsm, kf_cfg);
        init_estimates_none = get_initial_estimates(general_config_none, kf_cfg);
        
        [rx_sig_csm, los_phase, psi_csm] = get_received_signal(L1_C_over_N0_dBHz, 'CSM', doppler_profile, ...
            'S4', S4, 'tau0', tau0, 'simulation_time', simulation_time, 'settling_time', settling_time);
        [rx_sig_tppsm, ~, psi_tppsm] = get_received_signal(L1_C_over_N0_dBHz, 'TPPSM', doppler_profile, ...
            'tppsm_scenario', 'Severe', 'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', true);
        
        % Obtain state estimates for CSM
        [kf_ar_csm, ~]    = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_csm, 'CSM', adaptive_config_KF_AR);
        [akf_ar_csm, ~]   = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_csm, 'CSM', adaptive_config_AKF_AR);
        [ahl_kf_ar_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_csm, 'CSM', adaptive_config_AHL_KF_AR);
        [kf_std_csm, ~]   = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_none, 'none', adaptive_config_KF_std);
        [akf_std_csm, ~]  = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_none, 'none', adaptive_config_AKF_std);
        [ahl_kf_std_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_none, 'none', adaptive_config_AHL_KF_std);
        
        % Obtain state estimates for TPPSM
        [kf_ar_tppsm, ~]    = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_tppsm, 'TPPSM', adaptive_config_KF_AR);
        [akf_ar_tppsm, ~]   = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_tppsm, 'TPPSM', adaptive_config_AKF_AR);
        [ahl_kf_ar_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_tppsm, 'TPPSM', adaptive_config_AHL_KF_AR);
        [kf_std_tppsm, ~]   = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_none, 'none', adaptive_config_KF_std);
        [akf_std_tppsm, ~]  = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_none, 'none', adaptive_config_AKF_std);
        [ahl_kf_std_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_none, 'none', adaptive_config_AHL_KF_std);
        
        % Define validation vectors (assuming valid_time_vector is defined)
        valid_los_phase = wrapToPi(los_phase(valid_time_vector,1)); % This is being iterated unecessarily.
        valid_psi_phase_csm = angle(psi_csm(valid_time_vector,1));
        valid_psi_phase_tppsm = angle(psi_tppsm(valid_time_vector,1));
        valid_joint_phase_csm = wrapToPi(valid_los_phase + valid_psi_phase_csm);
        valid_joint_phase_tppsm = wrapToPi(valid_los_phase + valid_psi_phase_tppsm);
        
        % Store RMSE metrics in the results structure.
        % LOS Errors:
        results.csm.los.kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(kf_ar_csm(valid_time_vector,1)), valid_los_phase));
        results.csm.los.akf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(akf_ar_csm(valid_time_vector,1)), valid_los_phase));
        results.csm.los.ahl_kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(ahl_kf_ar_csm(valid_time_vector,1)), valid_los_phase));
        results.tppsm.los.kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(kf_ar_tppsm(valid_time_vector,1)), valid_los_phase));
        results.tppsm.los.akf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(akf_ar_tppsm(valid_time_vector,1)), valid_los_phase));
        results.tppsm.los.ahl_kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(ahl_kf_ar_tppsm(valid_time_vector,1)), valid_los_phase));
        
        % Scintillation Errors (only for AR estimates):
        results.csm.scint.kf_ar(m,pnv_idx) = sqrt(immse(kf_ar_csm(valid_time_vector, scint_index), valid_psi_phase_csm));
        results.csm.scint.akf_ar(m,pnv_idx) = sqrt(immse(akf_ar_csm(valid_time_vector, scint_index), valid_psi_phase_csm));
        results.csm.scint.ahl_kf_ar(m,pnv_idx) = sqrt(immse(ahl_kf_ar_csm(valid_time_vector, scint_index), valid_psi_phase_csm));
        results.tppsm.scint.kf_ar(m,pnv_idx) = sqrt(immse(kf_ar_tppsm(valid_time_vector, scint_index), valid_psi_phase_tppsm));
        results.tppsm.scint.akf_ar(m,pnv_idx) = sqrt(immse(akf_ar_tppsm(valid_time_vector, scint_index), valid_psi_phase_tppsm));
        results.tppsm.scint.ahl_kf_ar(m,pnv_idx) = sqrt(immse(ahl_kf_ar_tppsm(valid_time_vector, scint_index), valid_psi_phase_tppsm));
        
        % Joint Errors:
        results.csm.joint.kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(kf_ar_csm(valid_time_vector,1) + kf_ar_csm(valid_time_vector, scint_index)), valid_joint_phase_csm));
        results.csm.joint.akf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(akf_ar_csm(valid_time_vector,1) + akf_ar_csm(valid_time_vector, scint_index)), valid_joint_phase_csm));
        results.csm.joint.ahl_kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(ahl_kf_ar_csm(valid_time_vector,1) + ahl_kf_ar_csm(valid_time_vector, scint_index)),valid_joint_phase_csm));
        results.csm.joint.kf_std(m,pnv_idx) = sqrt(immse(wrapToPi(kf_std_csm(valid_time_vector,1)),valid_joint_phase_csm));
        results.csm.joint.akf_std(m,pnv_idx) = sqrt(immse(wrapToPi(akf_std_csm(valid_time_vector,1)), valid_joint_phase_csm));
        results.csm.joint.ahl_kf_std(m,pnv_idx) = sqrt(immse(wrapToPi(ahl_kf_std_csm(valid_time_vector,1)), valid_joint_phase_csm));
        
        results.tppsm.joint.kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(kf_ar_tppsm(valid_time_vector,1) + kf_ar_tppsm(valid_time_vector, scint_index)), valid_joint_phase_tppsm));
        results.tppsm.joint.akf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(akf_ar_tppsm(valid_time_vector,1) + akf_ar_tppsm(valid_time_vector, scint_index)), valid_joint_phase_tppsm));
        results.tppsm.joint.ahl_kf_ar(m,pnv_idx) = sqrt(immse(wrapToPi(ahl_kf_ar_tppsm(valid_time_vector,1) + ahl_kf_ar_tppsm(valid_time_vector, scint_index)), valid_joint_phase_tppsm));
        results.tppsm.joint.kf_std(m,pnv_idx) = sqrt(immse(wrapToPi(kf_std_tppsm(valid_time_vector,1)), valid_joint_phase_tppsm));
        results.tppsm.joint.akf_std(m,pnv_idx) = sqrt(immse(wrapToPi(akf_std_tppsm(valid_time_vector,1)), valid_joint_phase_tppsm));
        results.tppsm.joint.ahl_kf_std(m,pnv_idx) = sqrt(immse(wrapToPi(ahl_kf_std_tppsm(valid_time_vector,1)), valid_joint_phase_tppsm));
    end
end

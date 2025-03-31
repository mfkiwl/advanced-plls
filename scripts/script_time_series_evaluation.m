% Script: script_KFAR_estimates_comparison.m
%
% Description:
%   This script generates received signal data for both CSM and TPPSM scenarios
%   under severe scintillation, configures the Kalman PLL (KFAR) and standard KF 
%   algorithms with different adaptive update options, computes the 
%   corresponding state estimates, and then plots the comparison.
%
%   The following variants are computed:
%     - KF-AR      : KF with AR augmentation (no adaptive update)
%     - AKF-AR     : KF with AR augmentation using simplified adaptation (hard_limited = false)
%     - AHL-KF-AR  : KF with AR augmentation using simplified adaptation (hard_limited = true)
%     - KF-std     : Standard KF (no AR augmentation; training_scint_model = 'none', no adaptive update)
%     - AKF-std    : Standard KF using simplified adaptation (hard_limited = false)
%     - AHL-KF-std : Standard KF using simplified adaptation (hard_limited = true)
%
% [1] R. A. M. Lopes, F. Antreich, F. Fohlmeister, M. Kriegel and H. K. Kuga, "Ionospheric 
%     Scintillation Mitigation With Kalman PLLs Employing Radial Basis Function Networks," 
%     in IEEE Transactions on Aerospace and Electronic Systems, vol. 59, no. 5, pp. 6878-6893,
%     Oct. 2023, doi: 10.1109/TAES.2023.3281431
%  Author: Rodrigo de Lima Florindo
%  ORCID: https://orcid.org/0000-0003-0412-5583
%  Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', 'libs')));

% Main seed for generating the received signal and the training data set.
seed = 4;
rng(seed);

%% Generating the received signal for CSM and TPPSM under severe scintillation scenarios
doppler_profile = [0, 1000, 0.94];
sampling_interval = 0.01; % 100 Hz
L1_C_over_N0_dBHz = 42;
simulation_time = 300;
S4 = 0.8;
tau0 = 0.5;
settling_time = sampling_interval;
is_refractive_effects_removed_received_signal = false;
is_refractive_effects_removed_training_data = true;
is_unwrapping_used = false;
[rx_sig_csm, los_phase, psi_csm] = get_received_signal(L1_C_over_N0_dBHz, 'CSM', doppler_profile, ...
    'S4', S4, 'tau0', tau0, 'simulation_time', simulation_time, 'settling_time', settling_time);
[rx_sig_tppsm, ~, psi_tppsm, diffractive_phase_tppsm, refractive_phase_settled] = get_received_signal(L1_C_over_N0_dBHz, 'TPPSM', doppler_profile, ...
    'tppsm_scenario', 'Severe', 'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', is_refractive_effects_removed_received_signal);

%% Generating KF-AR configurations and obtaining initial estimates
cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');
  
training_simulation_time = 300;
training_data_config_csm = struct('scintillation_model', 'CSM', 'S4', S4, 'tau0', tau0, ...
                                  'simulation_time', training_simulation_time, 'sampling_interval', sampling_interval, 'is_unwrapping_used', is_unwrapping_used);
training_data_config_tppsm = struct('scintillation_model', 'TPPSM', 'scenario', 'Severe', ...
                                    'simulation_time', training_simulation_time, 'is_refractive_effects_removed', is_refractive_effects_removed_training_data, 'sampling_interval', sampling_interval, 'is_unwrapping_used', is_unwrapping_used);
training_data_config_none = struct('scintillation_model', 'none', 'sampling_interval', sampling_interval);

% Here, we used the same noise variance as used in [1, Section V; Subsection A]
process_noise_variance = 2.6*1e-8; 
ar_model_order = 5;
general_config_csm = struct( ...
  'discrete_wiener_model_config', { {1, 3, 0.01, [0, 0, process_noise_variance], 1} }, ...
  'scintillation_training_data_config', training_data_config_csm, ...
  'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
  'initial_states_distributions_boundaries', { {[-pi, pi], [-5, 5], [-0.1, 0.1]} }, ...
  'real_doppler_profile', doppler_profile, ...
  'augmentation_model_initializer', struct('id', 'aryule', 'model_params', struct('model_order', ar_model_order)), ...
  'is_use_cached_settings', false, ...
  'is_generate_random_initial_estimates', true, ...
  'is_enable_cmd_print', false ...
);

general_config_tppsm = general_config_csm;
general_config_tppsm.scintillation_training_data_config = training_data_config_tppsm;

general_config_none = general_config_csm;
general_config_none.scintillation_training_data_config = training_data_config_none;
general_config_none.augmentation_model_initializer.id = 'none';
general_config_none.augmentation_model_initializer.model_params = struct();
is_enable_cmd_print = true;

[~, init_estimates_csm] = get_kalman_pll_config(general_config_csm, cache_dir, is_enable_cmd_print);
[~, init_estimates_tppsm] = get_kalman_pll_config(general_config_tppsm, cache_dir, is_enable_cmd_print);
[kf_cfg, init_estimates_none] = get_kalman_pll_config(general_config_none, cache_dir, is_enable_cmd_print);

%% Define adaptive configuration structures
% For the simplified adaptive update, we require L1_C_over_N0_dBHz, sampling_interval, and threshold.
hard_limited_threshold = 38;  % Example threshold in dB; adjust as needed.
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

%% Online model learning configuration
% Rodrigo's Heuristics: 
% 1) It seems that using a settling time higher than sampling_interval makes it
% difficult to the KF-AR to use the online modules properly.
% 2) For a `window_size` less than 400, the AR phase estimates might be so 
% close to zero that the aryule raises an error saying that the window 
% data values must be finite.
% 3) Higher values of the window size (> 1500) makes the scintillation
% phase estimates looking more to the wrapped scintillation phase of the
% TPPSM.
% 4) Training the AR model with a unwrapped phase seems to not be a
% reasonable approach. The phase estimates in this case are nearly zero.
% 5) Using a very low value for the `process_noise_variance` makes it easier
% for the filters to diverge (test it for 2.6e-14, for an example). However,
% for higher process noise variance, the scintillation phase estimates
% becomes closer to zero.
online_mdl_learning_cfg = struct('is_online', false, 'learning_method', 'sliding_window', 'window_size', 1500);

%% Obtain state estimates for CSM
% For CSM, the training_scint_model is 'CSM' (AR augmented).
% 1. KF-AR     : No adaptive update.
% 2. AKF-AR    : Simplified adaptive update with hard_limited = false.
% 3. AHL-KF-AR : Simplified adaptive update with hard_limited = true.
[kf_ar_csm, ~]      = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_csm, 'CSM', adaptive_config_KF_AR, online_mdl_learning_cfg);
[akf_ar_csm, ~]  = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_csm, 'CSM', adaptive_config_AKF_AR, online_mdl_learning_cfg);
[ahl_kf_ar_csm, ~]   = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_csm, 'CSM', adaptive_config_AHL_KF_AR, online_mdl_learning_cfg);

% For standard KF estimates, training_scint_model is 'none'.
% 4. KF-std    : No adaptive update.
% 5. AKF-std   : Simplified adaptive update with hard_limited = false.
% 6. AHL-KF-std: Simplified adaptive update with hard_limited = true.
[kf_std_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_none, 'none', adaptive_config_KF_std, online_mdl_learning_cfg);
[akf_std_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_none, 'none', adaptive_config_AKF_std, online_mdl_learning_cfg);
[ahl_kf_std_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_none, 'none', adaptive_config_AHL_KF_std, online_mdl_learning_cfg);

%% Obtain state estimates for TPPSM
% For TPPSM, the training_scint_model is 'TPPSM' (AR augmented).
[kf_ar_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_tppsm, 'TPPSM', adaptive_config_KF_AR, online_mdl_learning_cfg);
[akf_ar_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_tppsm, 'TPPSM', adaptive_config_AKF_AR, online_mdl_learning_cfg);
[ahl_kf_ar_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_tppsm, 'TPPSM', adaptive_config_AHL_KF_AR, online_mdl_learning_cfg);

% For standard KF estimates with TPPSM received signal, training_scint_model is 'none'.
[kf_std_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_none, 'none', adaptive_config_KF_std, online_mdl_learning_cfg);
[akf_std_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_none, 'none', adaptive_config_AKF_std, online_mdl_learning_cfg);
[ahl_kf_std_tppsm, ~]  = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_none, 'none', adaptive_config_AHL_KF_std, online_mdl_learning_cfg);

%% Time vector generation
time_vector = sampling_interval:sampling_interval:simulation_time;

%% Plot the estimates comparison
is_save_figures = false;
%The plotting function should be modified to handle these inputs accordingly.
plot_time_series_full(...
    time_vector, ...
    los_phase, ...
    psi_csm, ...
    kf_ar_csm, ...      % KF-AR (no adaptive update, AR)
    akf_ar_csm, ...     % AKF-AR (simplified, hard_limited false, AR)
    ahl_kf_ar_csm, ...  % AHL-KF-AR (simplified, hard_limited true, AR)
    kf_std_csm, ...     % KF-std (no adaptive update, standard KF)
    akf_std_csm, ...    % AKF-std (simplified, hard_limited false, standard KF)
    ahl_kf_std_csm, ... % AHL-KF-std (simplified, hard_limited true, standard KF)
    kf_ar_tppsm, ...
    akf_ar_tppsm, ...
    ahl_kf_ar_tppsm, ...
    kf_std_tppsm, ...
    akf_std_tppsm, ...
    ahl_kf_std_tppsm, ...
    doppler_profile, ...
    L1_C_over_N0_dBHz, ...
    seed, ...
    process_noise_variance, ...
    ar_model_order, ...
    refractive_phase_settled, ...
    diffractive_phase_tppsm, ...
    is_refractive_effects_removed_received_signal, ...
    is_refractive_effects_removed_training_data, ...
    is_save_figures);

% From 100s to 150s
% zoomed_samples = 100/sampling_interval + 1 : 150/sampling_interval;
% zoomed_time_vector = time_vector(zoomed_samples);
% plot_time_series_zoom(...
%     zoomed_time_vector, ...
%     los_phase(zoomed_samples,:), ...
%     psi_csm(zoomed_samples,:), ...
%     kf_ar_csm(zoomed_samples,:), ...      % KF-AR (no adaptive update, AR)
%     akf_ar_csm(zoomed_samples,:), ...     % AKF-AR (simplified, hard_limited false, AR)
%     ahl_kf_ar_csm(zoomed_samples,:), ...  % AHL-KF-AR (simplified, hard_limited true, AR)
%     kf_std_csm(zoomed_samples,:), ...     % KF-std (no adaptive update, standard KF)
%     akf_std_csm(zoomed_samples,:), ...    % AKF-std (simplified, hard_limited false, standard KF)
%     ahl_kf_std_csm(zoomed_samples,:), ... % AHL-KF-std (simplified, hard_limited true, standard KF)
%     psi_tppsm(zoomed_samples,:), ...
%     kf_ar_tppsm(zoomed_samples,:), ...
%     akf_ar_tppsm(zoomed_samples,:), ...
%     ahl_kf_ar_tppsm(zoomed_samples,:), ...
%     kf_std_tppsm(zoomed_samples,:), ...
%     akf_std_tppsm(zoomed_samples,:), ...
%     ahl_kf_std_tppsm(zoomed_samples,:), ...
%     doppler_profile, ...
%     L1_C_over_N0_dBHz, ...
%     seed, ...
%     process_noise_variance, ...
%     ar_model_order, ...
%     ps_realization, ...
%     is_save_figures);

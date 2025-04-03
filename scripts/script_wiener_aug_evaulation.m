
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
[rx_sig_csm, los_phase, psi_csm] = get_received_signal(L1_C_over_N0_dBHz, 'CSM', doppler_profile, ...
    'S4', S4, 'tau0', tau0, 'simulation_time', simulation_time, 'settling_time', settling_time);
[rx_sig_tppsm, ~, psi_tppsm, diffractive_phase_tppsm, refractive_phase_settled] = get_received_signal(L1_C_over_N0_dBHz, 'TPPSM', doppler_profile, ...
    'tppsm_scenario', 'Severe', 'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', is_refractive_effects_removed_received_signal);

%% Generating KF-AR configurations and obtaining initial estimates
cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');

training_simulation_time = 300;
training_data_config = struct('scintillation_model', 'none', 'sampling_interval', sampling_interval);

% Here, we used the same noise variance as used in [1, Section V; Subsection A]
process_noise_variance_los = 1e-7; 
aug_wiener_mdl_order = 3;
process_noise_variance_aug = 1e-2;
general_config = struct( ...
  'discrete_wiener_model_config', { {1, 3, 0.01, [0, 0, process_noise_variance_los], 1} }, ...
  'scintillation_training_data_config', training_data_config, ...
  'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
  'initial_states_distributions_boundaries', { {[-pi, pi], [-0.1, 0.1], [-0.01, 0.01]} }, ...
  'real_doppler_profile', doppler_profile, ...
  'augmentation_model_initializer', struct('id', 'kinematic', 'model_params', struct('wiener_mdl_order', aug_wiener_mdl_order, 'process_noise_variance', process_noise_variance_aug)), ...
  'is_use_cached_settings', false, ...
  'is_generate_random_initial_estimates', true, ...
  'is_enable_cmd_print', false ...
);

is_enable_cmd_print = true;

[kf_cfg, init_estimates] = get_kalman_pll_config(general_config, cache_dir, is_enable_cmd_print);

%% Define adaptive configuration structures
adaptive_config_KF_std = struct(...
    'measurement_cov_adapt_algorithm', 'none', ...
    'states_cov_adapt_algorithm', 'none', ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', false));

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
% for heigher process noise variance, the scintillation phase estimates
% becomes closer to zero.
online_mdl_learning_cfg = struct('is_online', false);

%% Obtain state estimates for CSM
[kf_kinematic_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates, 'none', adaptive_config_KF_std, online_mdl_learning_cfg);
[kf_kinematic_tppsm, ~] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates, 'none', adaptive_config_KF_std, online_mdl_learning_cfg);

time_vector = sampling_interval:sampling_interval:simulation_time;

figure;
subplot(2,1,1);
plot(time_vector, [kf_kinematic_csm(:,1) - los_phase,kf_kinematic_csm(:,4),kf_kinematic_csm(:,1) + kf_kinematic_csm(:,4) - los_phase, angle(psi_csm)]);
legend({'LOS Phase error ', 'Kinematic Phase', 'Joint Phase Error', 'Unwrapped True Scint Phase'});
subplot(2,1,2);
plot(time_vector, [kf_kinematic_tppsm(:,1) - los_phase, kf_kinematic_tppsm(:,4), kf_kinematic_tppsm(:,1) + kf_kinematic_tppsm(:,4) - los_phase, unwrap(angle(psi_tppsm)), refractive_phase_settled]);
legend({'LOS phase Error', 'Kinematic Phase', 'Joint Phase', 'Unwrapped Scint phase', 'Refractive Phase'});
% 
% figure;
% plot(time_vector, unwrap(angle(psi_tppsm)) - refractive_phase_settled);

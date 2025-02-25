clear all; clc;

addpath(genpath(fullfile(pwd,'..','libs')));

%% Generating the received signal for CSM and TPPSM under sever scintillatio scenarios
doppler_profile = [0,0,0];
L1_C_over_N0_dBHz = 45;
[csm_rx_sig, csm_los_phase, csm_psi] = get_received_signal(L1_C_over_N0_dBHz, 'CSM', doppler_profile, ...
    'S4', 0.8, 'tau0', 0.7, 'simulation_time', 300, 'settling_time', 50);
[tppsm_rx_sig, tppsm_los_phase, tppsm_psi] = get_received_signal(L1_C_over_N0_dBHz, 'TPPSM', doppler_profile, ...
    'tppsm_scenario', 'Severe', 'simulation_time', 300, 'settling_time', 50, 'is_refractive_effects_removed', true);

%% Generating KF-AR configurations and getting their estimates for CSM and TPPSM under severe scintillation scenarios
cache_dir = [fileparts(mfilename('fullpath')), '\cache'];

csm_training_data_config = struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0',0.7, 'simulation_time', 900, 'sampling_interval',0.01);
tppsm_training_data_config = struct('scintillation_model', 'TPPSM', 'scenario', 'Severe', 'simulation_time', 900, 'is_refractive_effects_removed', true, 'sampling_interval', 0.01);

general_config_csm = struct( ...
  'discrete_wiener_model_config', { {1,3,0.01,[0,0,1e-2],1} }, ...
  'scintillation_training_data_config', csm_training_data_config, ...
  'var_minimum_order', 3, ...
  'var_maximum_order', 3, ...
  'C_over_N0_array_dBHz', 45, ...
  'initial_states_distributions_boundaries',{ {[-pi,pi], [-5,5], [-0.1,0.1]} }, ...
  'real_doppler_profile', doppler_profile, ...
  'is_use_cached_settings', true, ...
  'is_generate_random_initial_estimates', false ...
);

general_config_tppsm = general_config_csm;
general_config_tppsm.scintillation_training_data_config = tppsm_training_data_config;

% Update the CSM settings in the cache file, if it is empty and get the initial estimates.
[~, csm_init_estimates]= get_kalman_pll_config(general_config_csm, cache_dir);
% Update the CSM settings in the cache file, if it is empty, generate
% the variable kf_cfg with the CSM and TPPSM configs and obtain the initial estimates.
[kf_cfg, tppsm_init_estimates] = get_kalman_pll_config(general_config_tppsm, cache_dir);

[csm_state_estimates, csm_error_covariance_estimates] = get_kalman_pll_estimates(csm_rx_sig, kf_cfg, csm_init_estimates, 'CSM');
[tppsm_state_estimates, tppsm_error_covariance_estimates] = get_kalman_pll_estimates(tppsm_rx_sig, kf_cfg, tppsm_init_estimates, 'TPPSM');

bp = 1;
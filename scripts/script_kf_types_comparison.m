% script_kf_types_comparison.m
%
% This script serves to evaluate the behavior of different KF types under
% ionospheric scintillation events, which are synthetically generated using
% the TPPSM.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', 'libs')));

% Main seed for generating the received signal and the training data set.
seed = 6;
rng(seed);

%% Generating the received signal for the TPPSM
doppler_profile = [0, 0, 0]; % Synthetic Doppler profile
sampling_interval = 0.01; % 100 Hz
L1_C_over_N0_dBHz = 40;
simulation_time = 300;
settling_time = sampling_interval;
is_refractive_effects_removed_received_signal = false;
[rx_sig_tppsm, los_phase, psi_tppsm, diffractive_phase_tppsm, refractive_phase_settled] = get_received_signal(L1_C_over_N0_dBHz, 'TPPSM', doppler_profile, ...
    'tppsm_scenario', 'strong', 'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', is_refractive_effects_removed_received_signal);

%% Generating KF-AR configurations and obtaining initial estimates
cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');
training_data_config = struct('scintillation_model', 'none', 'sampling_interval', sampling_interval);

process_noise_variance_los = 2.6*1e-1; 
% The key feature will be another part of this struct that will
% characterize the adopted KF type as one of the following: 
% {'standard', 'extended', 'unscented', 'cubature'};
general_config_standard = struct( ...
  'kf_type', 'standard', ...
  'discrete_wiener_model_config', { {1, 3, 0.01, [0, 0, process_noise_variance_los], 1} }, ...
  'scintillation_training_data_config', training_data_config, ...
  'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
  'initial_states_distributions_boundaries', { {[-pi, pi], [-1, 1], [-0.01, 0.01]} }, ...
  'real_doppler_profile', doppler_profile, ...
  'augmentation_model_initializer', struct('id', 'none', 'model_params', struct()), ...
  'is_use_cached_settings', false, ...
  'is_generate_random_initial_estimates', false, ...
  'is_enable_cmd_print', false ...
);

% All configuration of the extended one is the same as the standard, except
% by its type.
general_config_extended = general_config_standard;
general_config_extended.kf_type = 'extended';

general_config_unscented = general_config_standard;
general_config_unscented.kf_type = 'unscented';

is_enable_cmd_print = true;

[~, ~] = get_kalman_pll_config(general_config_standard, cache_dir, is_enable_cmd_print);
[~, ~] = get_kalman_pll_config(general_config_extended, cache_dir, is_enable_cmd_print);
[kf_cfg, init_estimates] = get_kalman_pll_config(general_config_unscented, cache_dir, is_enable_cmd_print);

adaptive_cfg_nonadaptive = struct(...
    'measurement_cov_adapt_algorithm', 'none', ...
    'states_cov_adapt_algorithm', 'none', ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', false));

online_mdl_learning_cfg = struct('is_online', false);

[kf_std, error_cov_std] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates, 'standard', 'none', adaptive_cfg_nonadaptive, online_mdl_learning_cfg);
[kf_ext, error_cov_ext] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates, 'extended', 'none', adaptive_cfg_nonadaptive, online_mdl_learning_cfg);
[kf_uns, error_cov_uns] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates, 'unscented', 'none', adaptive_cfg_nonadaptive, online_mdl_learning_cfg);

time_vector = sampling_interval:sampling_interval:simulation_time;
true_total_phase = los_phase + get_corrected_phase(psi_tppsm);
phase_error_kf_std = kf_std(:,1) - true_total_phase;
phase_error_kf_ext = kf_ext(:,1) - true_total_phase;
phase_error_kf_uns = kf_uns(:,1) - true_total_phase;

linewidth = 2;

figure;
hold on;
scatter(time_vector,phase_error_kf_std, 'Marker', '.', 'LineWidth', linewidth);
scatter(time_vector,phase_error_kf_ext, 'Marker', '.', 'LineWidth', linewidth);
scatter(time_vector,phase_error_kf_uns, 'Marker', '.', 'LineWidth', linewidth);
legend({'Standard', 'Extended'}, Location="best");
hold off;

RMSE_standard = rms(wrapToPi(phase_error_kf_std));
RMSE_extended = rms(wrapToPi(phase_error_kf_ext));
RMSE_unscented = rms(wrapToPi(phase_error_kf_uns));

% Print the results to the MATLAB Command Window.
fprintf('RMSE (standard): %.4f\n', RMSE_standard);
fprintf('RMSE (extended): %.4f\n', RMSE_extended);
fprintf('RMSE (unscented): %.4f\n', RMSE_unscented);
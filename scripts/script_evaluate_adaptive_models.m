% script_evaluate_adaptive_models
%
% Description:
%
%  Author: Rodrigo de Lima Florindo
%  ORCID: https://orcid.org/0000-0003-0412-5583
%  Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', 'libs')));

% Main seed for generating the received signal and the training data set.
seed = 6;
rng(seed);

%% Generating the received signal for the TPPSM
doppler_profile = [0, 1000, 0.94]; % Synthetic Doppler profile
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
general_config = struct( ...
  'discrete_wiener_model_config', { {1, 3, 0.01, [0, 0, process_noise_variance_los], 1} }, ...
  'scintillation_training_data_config', training_data_config, ...
  'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
  'initial_states_distributions_boundaries', { {[-pi, pi], [-25, 25], [-0.01, 0.01]} }, ...
  'real_doppler_profile', doppler_profile, ...
  'augmentation_model_initializer', struct('id', 'none', 'model_params', struct()), ...
  'is_use_cached_settings', false, ...
  'is_generate_random_initial_estimates', true, ...
  'is_enable_cmd_print', false ...
);

is_enable_cmd_print = true;

[kf_cfg, init_estimates] = get_kalman_pll_config(general_config, cache_dir, is_enable_cmd_print);

%% Adaptive configurations to be evaluated
% NOTE: Hard limited is not considered for now.
hard_limited_flag = false;

% NWPR parameters
T_bit = 1/50;
M_nwpr = T_bit / sampling_interval;
N_nwpr = 20;

% Not adaptive: No measurement or state adaptation.
adaptive_cfg_nonadaptive = struct(...
    'measurement_cov_adapt_algorithm', 'none', ...
    'states_cov_adapt_algorithm', 'none', ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', hard_limited_flag));

% Use the NWPR approach for measurement covariance adaptation.
adaptive_cfg_nwpr = struct(...
    'measurement_cov_adapt_algorithm', 'nwpr', ...
    'measurement_cov_adapt_algorithm_params', struct('N_nwpr', N_nwpr, 'M_nwpr', M_nwpr), ...
    'states_cov_adapt_algorithm', 'none', ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', hard_limited_flag));

% Adapt both the measurements and the states covariance matrices:
% For the states adaptation using the 'matching' method.
window_size_matching = 25;
adaptive_cfg_nwpr_matching =struct('measurement_cov_adapt_algorithm', 'nwpr', ...
    'measurement_cov_adapt_algorithm_params', struct('N_nwpr', N_nwpr, 'M_nwpr', M_nwpr), ...
    'states_cov_adapt_algorithm', 'matching', ...
    'states_cov_adapt_algorithm_params', struct('method', 'RAE', 'window_size', window_size_matching), ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', false));

%% Online model learning configuration
% NOTE: Not applicable in this case, given that no augmentation model is
% being used. Thus, its `ìs_online` flag is configured as false.

online_mdl_learning_cfg = struct('is_online', false);

[kf, error_cov] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates, 'none', adaptive_cfg_nonadaptive, online_mdl_learning_cfg);
[kf_nwpr, error_cov_nwpr] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates, 'none', adaptive_cfg_nwpr, online_mdl_learning_cfg);
[kf_nwpr_matching, error_cov_nwpr_matching] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates, 'none', adaptive_cfg_nwpr_matching, online_mdl_learning_cfg);

time_vector = sampling_interval:sampling_interval:simulation_time;
true_total_phase = los_phase + get_corrected_phase(psi_tppsm);
phase_error_kf = kf(:,1) - true_total_phase;
phase_error_kf_nwpr = kf_nwpr(:,1) - true_total_phase;
phase_error_kf_nwpr_matching = kf_nwpr_matching(:,1) - true_total_phase;

figure;
hold on;
plot(time_vector,phase_error_kf, 'LineStyle', '-', 'LineWidth', 2);
plot(time_vector,phase_error_kf_nwpr, 'LineStyle', '-', 'LineWidth', 2);
plot(time_vector,phase_error_kf_nwpr_matching, 'LineStyle','-', 'LineWidth', 2);
legend({'KF', 'KF-NWPR', 'KF-NWPR-Matching'}, Location="best");
hold off;

time_vector = sampling_interval:sampling_interval:simulation_time;
diffractive_phase_error_kf = kf(:,1) - los_phase - refractive_phase_settled;
diffractive_phase_error_kf_nwpr = kf_nwpr(:,1) - los_phase - refractive_phase_settled;
diffractive_phase_error_kf_nwpr_matching = kf_nwpr_matching(:,1) - los_phase - refractive_phase_settled;

figure;
hold on;
plot(time_vector,diffractive_phase_error_kf, 'LineStyle', '-', 'LineWidth', 2);
plot(time_vector,diffractive_phase_error_kf_nwpr, 'LineStyle', '-', 'LineWidth', 2);
plot(time_vector,diffractive_phase_error_kf_nwpr_matching, 'LineStyle','-', 'LineWidth', 2);
plot(time_vector,get_corrected_phase(psi_tppsm) - refractive_phase_settled, 'LineStyle','-', 'LineWidth', 2, 'Color', 'k');
legend({'KF', 'KF-NWPR', 'KF-NWPR-Matching', 'Diffractive Phase'}, Location="best");
hold off;
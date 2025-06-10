% Script: script_arima_aug_evaluation.m
%
% Description:
%   This script simulates and compares the performance of Kalman PLL configurations
%   for mitigating ionospheric scintillation effects in two scenarios: CSM and TPPSM.
%   It generates the received signal under severe scintillation conditions, obtains
%   initial estimates and Kalman filter configurations (KF-ARIMA with ARIMA augmentation)
%   using predefined training and noise parameters, and then computes the state estimates.
%   Finally, the script visualizes the phase error and augmented phase estimates for both
%   scenarios.
%
%   The following main steps are performed:
%     1. Set up the environment and add required library paths.
%     2. Define simulation parameters such as seed, doppler profile, sampling interval,
%        C/N0, simulation time, and scintillation parameters.
%     3. Generate received signal data for the CSM and TPPSM scenarios using the function
%        get_received_signal.
%     4. Configure the Kalman PLL (KF-ARIMA) settings and obtain initial state estimates via
%        get_kalman_pll_config.
%     5. Set up adaptive configurations and online model learning parameters.
%     6. Compute state estimates for both scenarios using get_kalman_pll_estimates.
%     7. Plot the comparison of the phase errors, ARIMA phase, joint phase errors, and
%        true/unwrapped scintillation phases.
%
%  Author: Rodrigo de Lima Florindo
%  ORCID: https://orcid.org/0000-0003-0412-5583
%  Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', '..', 'libs')));

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
    'tppsm_scenario', 'strong', 'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', is_refractive_effects_removed_received_signal);

%% Generating KF-AR configurations and obtaining initial estimates
cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');

is_refractive_effects_removed_training_data = false;
is_unwrapping_used = true;
training_simulation_time = 300;
training_data_config_csm = struct('scintillation_model', 'CSM', ...
                                  'S4', S4, ...
                                  'tau0', tau0, ...
                                  'simulation_time', training_simulation_time, ...
                                  'sampling_interval', sampling_interval, ...
                                  'is_unwrapping_used', is_unwrapping_used);
training_data_config_tppsm = struct('scintillation_model', 'TPPSM', ...
                                    'scenario', 'strong', ...
                                    'simulation_time', training_simulation_time, ...
                                    'is_refractive_effects_removed', is_refractive_effects_removed_training_data, ...
                                    'sampling_interval', sampling_interval, ...
                                    'is_unwrapping_used', is_unwrapping_used);
% Here, we used the same noise variance as used in [1, Section V; Subsection A]
process_noise_variance_los = 1e-1; 
arima_p = 3;
arima_D = 1;
arima_q = 2;
general_config_csm = struct( ...
  'kf_type', 'standard', ...
  'discrete_wiener_model_config', { {1, 3, 0.01, [0, 0, process_noise_variance_los], 1} }, ...
  'scintillation_training_data_config', training_data_config_csm, ...
  'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
  'initial_states_distributions_boundaries', { {[-pi, pi], [-0.1, 0.1], [-0.01, 0.01]} }, ...
  'real_doppler_profile', doppler_profile, ...
  'augmentation_model_initializer', struct('id', 'arima', 'model_params', struct('p',arima_p,'D',arima_D,'q',arima_q)), ...
  'is_use_cached_settings', false, ...
  'is_generate_random_initial_estimates', false, ...
  'is_enable_cmd_print', false ...
);
general_config_tppsm = general_config_csm;
general_config_tppsm.scintillation_training_data_config = training_data_config_tppsm;

is_enable_cmd_print = true;

[~, init_estimates_csm] = get_kalman_pll_config(general_config_csm, cache_dir, is_enable_cmd_print);
[kf_cfg, init_estimates_tppsm] = get_kalman_pll_config(general_config_tppsm, cache_dir, is_enable_cmd_print);

%% Define adaptive configuration structures

T_bit = 1/50;
% Amount of samples within a navigation symbol period
M_nwpr = T_bit / sampling_interval;
% Amount of samples to be used to get an average C/N0 estimate
N_nwpr = 20;
adaptive_cfg = struct('measurement_cov_adapt_algorithm', 'nwpr', ...
                      'measurement_cov_adapt_algorithm_params', struct('N_nwpr', N_nwpr, 'M_nwpr', M_nwpr), ...
                      'states_cov_adapt_algorithm', 'none', ...
                      'sampling_interval', sampling_interval, ...
                      'hard_limited', struct('is_used', false));
% adaptive_cfg = struct('measurement_cov_adapt_algorithm', 'simplified', ...
%                       'measurement_cov_adapt_algorithm_params', struct('L1_C_over_N0_dBHz', L1_C_over_N0_dBHz), ...
%                       'states_cov_adapt_algorithm', 'none', ...
%                       'sampling_interval', sampling_interval, ...
%                       'hard_limited', struct('is_used', false));
%% Online model learning configuration

online_mdl_learning_cfg = struct('is_online', false);

%% Obtain state estimates for CSM
[kf_arima_csm, error_covariance_csm] = get_kalman_pll_estimates(rx_sig_csm, kf_cfg, init_estimates_csm, 'CSM', adaptive_cfg, online_mdl_learning_cfg);
[kf_arima_tppsm, error_covariance_tppsm, L1_c_over_n0_linear_estimates_tppsm] = get_kalman_pll_estimates(rx_sig_tppsm, kf_cfg, init_estimates_tppsm, 'TPPSM', adaptive_cfg, online_mdl_learning_cfg);

time_vector = sampling_interval:sampling_interval:simulation_time;

arima_est_csm = kf_arima_csm(:,4) + kf_arima_csm(:,5);
arima_est_tppsm = kf_arima_tppsm(:,4) + kf_arima_tppsm(:,5);
figure;
subplot(2,1,1);
plot(time_vector, [kf_arima_csm(:,1) - los_phase, arima_est_csm,kf_arima_csm(:,1) + arima_est_csm - los_phase, unwrap(angle(psi_csm))]);
legend({'LOS Phase error ', 'ARIMA Phase', 'Joint Phase Error', 'Unwrapped True Scint Phase'});
subplot(2,1,2);
plot(time_vector, [kf_arima_tppsm(:,1) - los_phase, arima_est_tppsm, kf_arima_tppsm(:,1) + arima_est_tppsm - los_phase, unwrap(angle(psi_tppsm)), refractive_phase_settled]);
legend({'LOS Phase Error', 'ARIMA Phase', 'Joint Phase Error', 'Unwrapped Scint phase', 'Refractive Phase'});


trace_ts = zeros(size(error_covariance_tppsm,1),1);
for i = 1:size(error_covariance_tppsm,1)
    trace_ts(i) = trace(squeeze(error_covariance_tppsm(i,:,:)));
end

figure;

plot(time_vector, trace_ts);

L1_C_over_N0_dbhz_estimates_tppsm = 10*log10(L1_c_over_n0_linear_estimates_tppsm);
true_carrier_to_noise_ratio_dbhz = 10*log10(abs(rx_sig_tppsm).^2) + L1_C_over_N0_dBHz;
figure;
plot(time_vector, [true_carrier_to_noise_ratio_dbhz,L1_C_over_N0_dbhz_estimates_tppsm]);

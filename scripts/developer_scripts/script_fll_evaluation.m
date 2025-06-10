%% FLL Evaluation Script for Thermal Noise Only
% script_fll_evaluation.m
%
% [Description Placeholder]
%
%  Author: Rodrigo de Lima Florindo
%  ORCID: https://orcid.org/0000-0003-0412-5583
%  Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', '..', 'libs')));

% Set simulation parameters
doppler_profile = [0, 1000, 0.94];
sampling_interval = 0.01; % 100 Hz
L1_C_over_N0_dBHz = 42;
simulation_time = 300;

% Thermal noise only: no scintillation effects.
% los_phase is the true phase (in radians) due to Doppler/LOS dynamics.
[rx_sig, los_phase] = get_received_signal(L1_C_over_N0_dBHz, 'none', ...
    doppler_profile, 'simulation_time', simulation_time);

%% Define FLL configuration sets

% Configuration 1: First order loop filter
fll_config_1 = struct;
fll_config_1.discriminator_type = 'atan2';
fll_config_1.sampling_interval = sampling_interval;
fll_config_1.initial_frequency_estimate_boundaries = [-0.1, 0.1];
fll_config_1.loop_filter_cfg = struct('order', 1, 'noise_bandwidth', 10, 'kd', 1);

% Configuration 2: Second order loop filter (critically damped)
fll_config_2 = struct;
fll_config_2.discriminator_type = 'atan2';
fll_config_2.sampling_interval = sampling_interval;
fll_config_2.initial_frequency_estimate_boundaries = [-0.1, 0.1];
fll_config_2.loop_filter_cfg = struct('order', 2, 'noise_bandwidth', 10, ...
    'damping_mode', 'critically_damped', 'kd', 1);

% Configuration 3: Second order loop filter (underdamped)
fll_config_3 = struct;
fll_config_3.discriminator_type = 'atan2';
fll_config_3.sampling_interval = sampling_interval;
fll_config_3.initial_frequency_estimate_boundaries = [-0.1, 0.1];
fll_config_3.loop_filter_cfg = struct('order', 2, 'noise_bandwidth', 10, ...
    'damping_mode', 'underdamped', 'kd', 1);

%% Get FLL estimates for thermal noise (no scintillation)
[phase_est_1, freq_est_1] = get_fll_estimates(fll_config_1, rx_sig);
[phase_est_2, freq_est_2] = get_fll_estimates(fll_config_2, rx_sig);
[phase_est_3, freq_est_3] = get_fll_estimates(fll_config_3, rx_sig);

%% Plot the results: Compare FLL estimated phase with the true LOS phase
figure;
plot(los_phase, 'k', 'LineWidth', 1.5); hold on;
plot(phase_est_1, 'b--');
plot(phase_est_2, 'r-.');
plot(phase_est_3, 'g:');
hold off;
title('Thermal Noise: Phase Estimates vs. True LOS Phase');
xlabel('Sample');
ylabel('Phase (radians)');
legend('True LOS Phase', '1st Order', '2nd Order, Critically Damp.', ...
       '2nd Order, Underdamp.', 'Location', 'Best');

%% (Optional) Plot Frequency Estimates
figure;
plot(freq_est_1, 'b--'); hold on;
plot(freq_est_2, 'r-.');
plot(freq_est_3, 'g:');
hold off;
title('Thermal Noise: Frequency Estimates for Different FLL Configurations');
xlabel('Sample');
ylabel('Frequency (Hz)');
legend('1st Order', '2nd Order, Critically Damp.', '2nd Order, Underdamp.', 'Location', 'Best');

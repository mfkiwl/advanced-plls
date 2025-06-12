% Script: script_KFAR_ar_phase_estimates_comparison.m
%
% Description:
%
%
% [1] R. A. M. Lopes, F. Antreich, F. Fohlmeister, M. Kriegel and H. K. Kuga, "Ionospheric 
%     Scintillation Mitigation With Kalman PLLs Employing Radial Basis Function Networks," 
%     in IEEE Transactions on Aerospace and Electronic Systems, vol. 59, no. 5, pp. 6878-6893,
%     Oct. 2023, doi: 10.1109/TAES.2023.3281431
%  Author: Rodrigo de Lima Florindo
%  ORCID: https://orcid.org/0000-0003-0412-5583
%  Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', '..', 'libs')));

% Main seed for generating the received signal and the training data set.
seed = 3;
rng(seed);

%% Generating the received signal for CSM and CPSSM under scintillation scenarios
doppler_profile = [0, 1000, 0.94];
sampling_interval = 0.01; % 100 Hz
L1_C_over_N0_dBHz = 42;
simulation_time = 300;
settling_time = 50;
cpssm_scenario = 'weak';
is_refractive_effects_removed_received_signal = false; % Change this to false to assess the behavior with refractive effects, and to true to assess only without refractive effects.
is_refractive_effects_removed_training_data = true;
is_unwrapping_used = false;
[rx_sig_cpssm, ~, psi_cpssm, diffractive_phase_cpssm, refractive_phase_settled] = get_received_signal(L1_C_over_N0_dBHz, 'TPPSM', doppler_profile, seed, ...
    'tppsm_scenario', cpssm_scenario, 'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', is_refractive_effects_removed_received_signal);

%% Generating KF-AR configurations and obtaining initial estimates
cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');
  
training_simulation_time = 300;
training_data_config_cpssm = struct( ...
    'scintillation_model', 'TPPSM', ...
    'scenario', cpssm_scenario, ...
    'simulation_time', training_simulation_time, 'is_refractive_effects_removed', ...
    is_refractive_effects_removed_training_data, 'sampling_interval', ...
    sampling_interval, 'is_unwrapping_used', is_unwrapping_used);


sigma2_W_3_amount = 5;
sigma2_W_3_sweep = logspace(-8,2, sigma2_W_3_amount);

ar_model_order = 14;
expected_doppler_profile = [0,1000,0.94];

ar_phase_time_series = zeros((simulation_time - settling_time)/sampling_interval, sigma2_W_3_amount);
is_enable_cmd_print = true;

%% Define configuration structure
adaptive_config_KF_AR = struct(...
    'measurement_cov_adapt_algorithm', 'none', ...
    'states_cov_adapt_algorithm', 'none', ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', false));

%% Online model learning configuration
online_mdl_learning_cfg = struct('is_online', false);

for sigma2_W_3_idx = 1:sigma2_W_3_amount
    general_config_cpssm = struct( ...
        'kf_type', 'standard', ...
        'discrete_wiener_model_config', { {1, 3, 0.01, [0, 0, sigma2_W_3_sweep(sigma2_W_3_idx)], 1} }, ...
        'scintillation_training_data_config', training_data_config_cpssm, ...
        'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
        'initial_states_distributions_boundaries', { {[-pi, pi], [-25, 25], [-1, 1]} }, ...
        'expected_doppler_profile', expected_doppler_profile, ...
        'augmentation_model_initializer', struct('id', 'aryule', 'model_params', struct('model_order', ar_model_order)), ...
        'is_use_cached_settings', false, ...
        'is_generate_random_initial_estimates', true, ...
        'is_enable_cmd_print', false ...
    );
    [kf_cfg, init_estimates_cpssm] = get_kalman_pll_config(general_config_cpssm, cache_dir, is_enable_cmd_print);
    [kf_ar_cpssm, ~] = get_kalman_pll_estimates(rx_sig_cpssm, kf_cfg, init_estimates_cpssm, 'standard', 'TPPSM', adaptive_config_KF_AR, online_mdl_learning_cfg);
    ar_phase_time_series(:,sigma2_W_3_idx) = kf_ar_cpssm(settling_time/sampling_interval + 1:end, 4); % L = 1 (Amount of bands) M = 3 (Wiener order); L + M -> AR phase estimates
end

% Zoomed-in AR vs. CPSSM phase plot (100–120 s)

% 1) Define sample & time windows
start_t = 200;      % seconds
end_t   = 220;      % seconds
% compute integer sample indices for full-rate signal
idx_full = round(start_t/sampling_interval) : round(end_t/sampling_interval);
% matching time vector
time_vector = (start_t : sampling_interval : end_t).';

% 2) Create figure with editable size (cm)
figWidth  = 40;
figHeight = 20;
fig = figure('Units','centimeters', ...
             'Position',[5 5 figWidth figHeight]);

hold on;

% 3) Plot CPSSM (diffractive) phase
h0 = plot( ...
    time_vector, ...
    diffractive_phase_cpssm(idx_full), ...
    'LineWidth', 2.5, ...
    'Color', 'k', ...
    'DisplayName', 'Diffractive Phase CPSSM' ...
);

% 4) Prepare 'cool' colormap for AR traces
colors = cool(sigma2_W_3_amount);

% 5) Plot each KF-AR series (note AR time series starts at sample 5000)
for idx = 1:sigma2_W_3_amount
    ar_idx = idx_full - 5000 + 1;  % map full-signal samples to AR-series indices
    plot( ...
        time_vector, ...
        ar_phase_time_series(ar_idx, idx), ...
        'LineWidth', 1.5, ...
        'Color', colors(idx,:), ...
        'DisplayName', sprintf('$\\sigma^2_{\\mathrm{W},3}=%.3g$', sigma2_W_3_sweep(idx)) ...
    );
end

hold off;

% 6) Axes labels (LaTeX + Times)
xlabel('Time [s]', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman');
ylabel('KF-AR - $\hat{\phi}_{\mathrm{AR},1}[k \mid k-1]$', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman');

% 7) Grids
grid on;
grid minor;

% 8) Legend
legend('Location','best', ...
       'Interpreter','latex', ...
       'FontName','Times New Roman');

% 9) Axis formatting
set(gca, ...
    'FontName','Times New Roman', ...
    'TickLabelInterpreter','latex', ...
    'FontSize', 22 ...
);

% 10) Ensure results folder exists & save
resultsDir = fullfile(pwd,'results');
if ~exist(resultsDir,'dir')
    mkdir(resultsDir);
end
exportgraphics(fig, fullfile(resultsDir,'ar_phase_zoom_100_120.pdf'), 'ContentType','vector');
savefig(fig,    fullfile(resultsDir,'ar_phase_zoom_100_120.fig'));

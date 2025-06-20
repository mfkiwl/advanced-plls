% plot_tppsm_strong_phase_error_kf.m
% Plots the total carrier phase error for the TPPSM under strong scintillation
% using only the standard Kalman Filter (KF) approach.

clearvars; clc;
addpath(genpath(fullfile(pwd, '..', '..', '..', 'libs')));

% 1) Simulation settings
seed               = 12;
rng(seed);
sampling_interval  = 1e-2;    % 100 Hz
simulation_time    = 300;     % seconds
settling_time      = 50;      % seconds
L1_C_over_N0_dBHz  = 42;      % dB-Hz
doppler_profile    = [0, 0, 0];

% 2) Generate TPPSM received signal under strong scintillation
is_refractive_removed = false;
[rx, true_los, psi_settled, ~, refractive_phase] = ...
    get_received_signal(L1_C_over_N0_dBHz, 'TPPSM', doppler_profile, seed, ...
                        'tppsm_scenario', 'strong', ...
                        'simulation_time', simulation_time, ...
                        'settling_time', settling_time, ...
                        'rhof_veff_ratio', 0.5, ...
                        'is_refractive_effects_removed', is_refractive_removed);
psi_phase = get_corrected_phase(psi_settled);
% 3) Configure standard KF (no AR augmentation)
sigma2_W_3 = 1e-5;
discrete_cfg = {1, 3, sampling_interval, [0, 0, sigma2_W_3], 1};
training_data_config_none = struct('scintillation_model', 'none', 'sampling_interval', sampling_interval);

gen_cfg = struct( ...
    'kf_type', 'standard', ...
    'discrete_wiener_model_config', {discrete_cfg}, ...
    'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
    'expected_doppler_profile', doppler_profile, ...
    'initial_states_distributions_boundaries', { {[-pi, pi], [-25, 25], [-1, 1]} }, ...
    'augmentation_model_initializer', struct('id', 'none', 'model_params', struct()), ...
    'scintillation_training_data_config', training_data_config_none, ...
    'is_use_cached_settings', false, ...
    'is_generate_random_initial_estimates', false, ...
    'is_enable_cmd_print', false ...
);

cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');
[gen_cfg, init_estimates] = get_kalman_pll_config(gen_cfg, cache_dir, false);

kf_cfg = struct( ...
    'measurement_cov_adapt_algorithm', 'none', ...
    'states_cov_adapt_algorithm', 'none', ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', false) ...
);
online_cfg = struct('is_online', false);

% 4) Compute KF estimates and error
[estimates, ~] = get_kalman_pll_estimates(rx, gen_cfg, init_estimates, ...
                    'standard', 'none', kf_cfg, online_cfg);
phi_w = estimates(:,1);
phi_t = phi_w;

% 4) Time vector and zoom indices
time = (0 : sampling_interval : simulation_time).';
start_t = 50; end_t = 250;
idx = time>=start_t & time<=end_t;

% 5) Plot 2×1 layout
figure('Color','w', 'Position', [50,50,1000,500]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% Tile 1: Estimates vs Ionospheric phase
nexttile; hold on;
plot(time(idx), phi_t(idx), 'b-', 'LineWidth',2, 'DisplayName','$\hat{\phi}_{\mathrm{W}, 1} [k \mid k - 1]$');
plot(time(idx), true_los(idx)+psi_phase(idx), 'g-','LineWidth',2,'DisplayName','$\phi_{\mathrm{I}, 1} [k]$');
plot(time(idx), refractive_phase(idx), 'k-','LineWidth',2,'DisplayName','$\phi_{\mathrm{R}, 1} [k]$');
ylabel('Phase [rad]','Interpreter','latex');
legend('Location','best','Interpreter','latex'); grid on; grid minor;
hold off;
set(gca, 'FontSize', 13);

% Tile 2: Refractive phase error
nexttile; hold on;
refr_error = phi_t(idx) - (refractive_phase(idx));
plot(time(idx), refr_error, 'b-','LineWidth',2,'DisplayName','$\hat{\phi}_{\mathrm{W}, 1} [k \mid k - 1] - \phi_{\mathrm{R}, 1} [k]$');
plot(time(idx), psi_phase(idx) - refractive_phase(idx), 'g-','LineWidth',2,'DisplayName','$\phi_{\mathrm{D}, 1} [k] = \phi_{\mathrm{I}, 1} [k] - \phi_{\mathrm{R}, 1} [k]$');
ylabel('Phase error [rad]','Interpreter','latex');
legend('Location','best','Interpreter','latex'); grid on; grid minor;
hold off;
set(gca, 'FontSize', 13);

% Tile 3: Total phase error
nexttile; hold on;
total_error = phi_t(idx) - (true_los(idx)+psi_phase(idx));
plot(time(idx), total_error, 'b-','LineWidth',2,'DisplayName','$\hat{\phi}_{\mathrm{W}, 1} [k \mid k - 1] - \phi_{\mathrm{I}, 1} [k]$');
ylabel('Phase error [rad]','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
legend('Location','best','Interpreter','latex'); grid on; grid minor;
hold off;
set(gca, 'FontSize', 13);

% 6) Save figure
out = fullfile('results', sprintf('cycle_slip_example'));
exportgraphics(gcf, [out,'.pdf'],'ContentType','vector');
savefig([out,'.fig']);

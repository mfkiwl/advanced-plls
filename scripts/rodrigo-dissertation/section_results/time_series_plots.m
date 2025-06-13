% script_plot_kfar_phases_zoom.m
%
% Template for plotting KF-AR phase components (φ_T, φ_W, φ_AR)
% in a 3×2 tile layout (rows=phase, cols=severity) for each model.
% Zoom window: 150–170 s, single realization (seed=3).
%
% Author: [Your Name]
% Date:   [YYYY-MM-DD]

clearvars; clc;

% add libraries
addpath(genpath(fullfile(pwd,'..','..','..','libs')));

%% Reproducibility
seed = 2;
rng(seed);

%% Configuration
cache_dir       = fullfile(fileparts(mfilename('fullpath')),'cache');
sampling_interval = 1e-2;
settling_time     = 50;
simulation_time   = 300;

% get only KF-AR config
[ahl_kf_ar_cfg, online_mdl_learning_cfg] = get_adaptive_cfgs();

%% Zoom window in seconds & samples
zoom_start = 50;    % seconds
zoom_end   = 100;    % seconds
idx_full   = round(zoom_start/sampling_interval) : round(zoom_end/sampling_interval);
time_zoom  = (zoom_start: sampling_interval : zoom_end).';

%% Models, Severities, Phases
model_names = {'csm', 'cpssm_wo_refr', 'cpssm_w_refr'};
severities  = {'weak', 'strong'};
field_names = {'phi_T','phi_W','phi_AR'};
field_labels = { ...
    '$\hat{\phi}_{\mathrm{T},1} [k \mid k - 1] - \phi_{\mathrm{T},1} [k]$', ...
    '$\hat{\phi}_{\mathrm{W},1} [k \mid k - 1] - \phi_{\mathrm{LOS},1} [k]$', ...
    '$\hat{\phi}_{\mathrm{AR},1} [k \mid k - 1] - \phi_{\mathrm{D},1} [k]$', ...
};

results_dir = fullfile(pwd,'results');
if ~exist(results_dir,'dir')
    mkdir(results_dir);
end

sigma2_W_3_amount = 5;
sigma2_W_3_sweep = logspace(-2,-10, sigma2_W_3_amount);

%% loop over models

% CSM
fig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9],'Color','w');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for row = 1:3
    for col = 1:2
        ax = nexttile;
        hold(ax,'on');

        severity = severities{col};

        % --- generate received signal & diffractive phase once ---
        [rx_inputs_base, ~, init_csm, ~, ~, ar_phase_idx] = ...
            get_overall_cfgs(cache_dir, 'csm', severity, true, sigma2_W_3_sweep(1), ...
                                sampling_interval, settling_time, simulation_time, seed);
        [rx_signal, true_los_phase, ~, diffractive_phase] = ...
            get_received_signal(rx_inputs_base{:});

        
        % prepare colors
        colors = cool(sigma2_W_3_amount);

        % loop over σ²_{W,3}
        for k = 1:sigma2_W_3_amount
            sigma2_W_3 = sigma2_W_3_sweep(k);

            % get KF-AR config & init estimates for this σ²
            [~, gen_kf_cfg, init_csm, init_cpssm, ~, ~] = ...
                get_overall_cfgs(cache_dir, 'csm', severity, true, sigma2_W_3, ...
                                    sampling_interval, settling_time, simulation_time, seed);

            % compute KF-AR estimates
            [kf_ar_est, ~] = get_kalman_pll_estimates( ...
                rx_signal, gen_kf_cfg, init_csm, ...
                'standard', upper('csm'), ahl_kf_ar_cfg, online_mdl_learning_cfg);

            switch row
                case 1  % total phase error
                    phi_T = kf_ar_est(:,1) + kf_ar_est(:,ar_phase_idx);
                    phase_error = phi_T(idx_full) - (true_los_phase(idx_full) + diffractive_phase(idx_full));
                case 2  % Wiener noise error
                    phi_W = kf_ar_est(:,1);
                    phase_error = phi_W(idx_full) - true_los_phase(idx_full);
                case 3  % AR component error
                    phi_AR = kf_ar_est(:,ar_phase_idx);
                    phase_error = phi_AR(idx_full) - diffractive_phase(idx_full);
            end

            % plot error for this σ²
            plot(ax, time_zoom, phase_error, ...
                'LineWidth',1.5, ...
                'Color',colors(k,:), ...
                'DisplayName',sprintf('$\\sigma^2_{W,3}=%.1e$',sigma2_W_3));
        end

        % plot diffractive_phase error = 0 baseline (in black)
        reference_phase_error = zeros(size(idx_full));

        % if row == 3
        %     reference_phase_error = diffractive_phase(idx_full);
        % end
        plot(ax, time_zoom, reference_phase_error, ...
            'LineWidth',2.5,'Color','k', ...
            'DisplayName', 'Target Error');

        hold(ax,'off');

        % only add legend on first row, first column
        if row == 1 && col == 1
            legend(ax, 'Location','best', ...
                'Interpreter','latex', ...
                'FontName','Times New Roman');
        end

        if col == 1
            severity_title = 'Weak';
        elseif col == 2
            severity_title = 'Strong';
        end

        % titles & labels
        if row==1
            title(ax, severity_title, 'Interpreter','latex','FontName','Times New Roman');
        end
        if col==1
            ylabel(ax, field_labels{row}, 'Interpreter','latex','FontName','Times New Roman');
        end
        if row==3
            xlabel(ax, 'Time [s]', 'Interpreter','latex','FontName','Times New Roman');
        end

        % styling
        grid(ax,'on'); grid(ax,'minor');
        set(ax,'FontName','Times New Roman','TickLabelInterpreter','latex', 'FontSize', 22);
    end
end

% save outputs
out_base = fullfile(results_dir, sprintf('%s_phase_errors_zoom','csm'));
exportgraphics(fig, [out_base,'.pdf'], 'ContentType','vector');
savefig(fig, [out_base,'.fig']);

% CPSSM without refractive phase
fig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9],'Color','w');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for row = 1:3
    for col = 1:2
        ax = nexttile;
        hold(ax,'on');

        severity = severities{col};

        % --- generate received signal & diffractive phase once ---
        [rx_inputs_base, ~, ~, init_cpssm, ~, ar_phase_idx] = ...
            get_overall_cfgs(cache_dir, 'cpssm', severity, true, sigma2_W_3_sweep(1), ...
                                sampling_interval, settling_time, simulation_time, seed);
        [rx_signal, true_los_phase, ~, diffractive_phase] = ...
            get_received_signal(rx_inputs_base{:});
        
        % prepare colors
        colors = cool(sigma2_W_3_amount);

        % loop over σ²_{W,3}
        for k = 1:sigma2_W_3_amount
            sigma2_W_3 = sigma2_W_3_sweep(k);

            % get KF-AR config & init estimates for this σ²
            [~, gen_kf_cfg, init_csm, init_cpssm, ~, ~] = ...
                get_overall_cfgs(cache_dir, 'cpssm', severity, true, sigma2_W_3, ...
                                    sampling_interval, settling_time, simulation_time, seed);

            % Compute KF-AR estimates
            [kf_ar_est, ~] = get_kalman_pll_estimates( ...
                rx_signal, gen_kf_cfg, init_csm, ...
                'standard', upper('tppsm'), ahl_kf_ar_cfg, online_mdl_learning_cfg);

            switch row
                case 1  % total phase error
                    phi_T = kf_ar_est(:,1) + kf_ar_est(:,ar_phase_idx);
                    phase_error = phi_T(idx_full) - (true_los_phase(idx_full) + diffractive_phase(idx_full));
                case 2  % Wiener noise error
                    phi_W = kf_ar_est(:,1);
                    phase_error = phi_W(idx_full) - true_los_phase(idx_full);
                case 3  % AR component error
                    phi_AR = kf_ar_est(:,ar_phase_idx);
                    phase_error = phi_AR(idx_full) - diffractive_phase(idx_full);
            end

            % plot error for this σ²
            plot(ax, time_zoom, phase_error, ...
                'LineWidth',1.5, ...
                'Color',colors(k,:), ...
                'DisplayName',sprintf('$\\sigma^2_{W,3}=%.1e$',sigma2_W_3));
        end

        % plot diffractive_phase error = 0 baseline (in black)
        reference_phase_error = zeros(size(idx_full));

        % if row == 3
        %     reference_phase_error = diffractive_phase(idx_full);
        % end
        plot(ax, time_zoom, reference_phase_error, ...
            'LineWidth',2.5,'Color','k', ...
            'DisplayName', 'Target Error');

        hold(ax,'off');

        % only add legend on first row, first column
        if row == 1 && col == 1
            legend(ax, 'Location','best', ...
                'Interpreter','latex', ...
                'FontName','Times New Roman');
        end

        if col == 1
            severity_title = 'Weak';
        elseif col == 2
            severity_title = 'Strong';
        end

        % titles & labels
        if row==1
            title(ax, severity_title, 'Interpreter','latex','FontName','Times New Roman');
        end
        if col==1
            ylabel(ax, field_labels{row}, 'Interpreter','latex','FontName','Times New Roman');
        end
        if row==3
            xlabel(ax, 'Time [s]', 'Interpreter','latex','FontName','Times New Roman');
        end

        % styling
        grid(ax,'on'); grid(ax,'minor');
        set(ax,'FontName','Times New Roman','TickLabelInterpreter','latex', 'FontSize', 22);
    end
end

% CPSSM with refractive phase component
fig = figure('Units','normalized','Position',[0.05 0.05 0.9 0.9],'Color','w');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for row = 1:3
    for col = 1:2
        ax = nexttile;
        hold(ax,'on');

        severity = severities{col};

        % --- generate received signal & diffractive phase once ---
        [rx_inputs_base, ~, ~, init_cpssm, ~, ar_phase_idx] = ...
            get_overall_cfgs(cache_dir, 'cpssm', severity, false, sigma2_W_3_sweep(1), ...
                                sampling_interval, settling_time, simulation_time, seed);
        [rx_signal, true_los_phase, psi_cpssm, diffractive_phase, refractive_phase] = ...
            get_received_signal(rx_inputs_base{:});
        
        % prepare colors
        colors = cool(sigma2_W_3_amount);

        % loop over σ²_{W,3}
        for k = 1:sigma2_W_3_amount
            sigma2_W_3 = sigma2_W_3_sweep(k);

            % get KF-AR config & init estimates for this σ²
            [~, gen_kf_cfg, init_csm, init_cpssm, ~, ~] = ...
                get_overall_cfgs(cache_dir, 'cpssm', severity, false, sigma2_W_3, ...
                                    sampling_interval, settling_time, simulation_time, seed);

            % Compute KF-AR estimates
            [kf_ar_est, ~] = get_kalman_pll_estimates( ...
                rx_signal, gen_kf_cfg, init_csm, ...
                'standard', upper('tppsm'), ahl_kf_ar_cfg, online_mdl_learning_cfg);

            full_phase_propagated_scint_field = get_corrected_phase(psi_cpssm);
            switch row
                case 1  % total phase error
                    phi_T = kf_ar_est(:,1) + kf_ar_est(:,ar_phase_idx);
                    phase_error = phi_T(idx_full) - (true_los_phase(idx_full) + full_phase_propagated_scint_field(idx_full));
                case 2  % Wiener noise error
                    phi_W = kf_ar_est(:,1);
                    phase_error = phi_W(idx_full) - true_los_phase(idx_full);
                case 3  % AR component error
                    phi_AR = kf_ar_est(:,ar_phase_idx);
                    phase_error = phi_AR(idx_full) - diffractive_phase(idx_full);
            end

            % plot error for this σ²
            plot(ax, time_zoom, phase_error, ...
                'LineWidth',1.5, ...
                'Color',colors(k,:), ...
                'DisplayName',sprintf('$\\sigma^2_{W,3}=%.1e$',sigma2_W_3));
        end

        % plot diffractive_phase error = 0 baseline (in black)
        reference_phase_error = zeros(size(idx_full));
        
        % if row == 3
        %     reference_phase_error = diffractive_phase(idx_full);
        % end
        plot(ax, time_zoom, reference_phase_error, ...
            'LineWidth',2.5,'Color','k', ...
            'DisplayName', 'Target Error');

        if row == 2
            plot(ax, time_zoom, refractive_phase(idx_full), ...
                'LineWidth', 2.5,'Color','k', ...
                'LineStyle', ':', ...
                'DisplayName', 'Refractive Phase ($\phi_{\mathrm{R}, 1} [k]$)');
        end

        hold(ax,'off');

        % only add legend on first row, first column
        if row == 2 && col == 1
            legend(ax, 'Location','best', ...
                'Interpreter','latex', ...
                'FontName','Times New Roman');
        end

        if col == 1
            severity_title = 'Weak';
        elseif col == 2
            severity_title = 'Strong';
        end

        % titles & labels
        if row==1
            title(ax, severity_title, 'Interpreter','latex','FontName','Times New Roman');
        end
        if col==1
            ylabel(ax, field_labels{row}, 'Interpreter','latex','FontName','Times New Roman');
        end
        if row==3
            xlabel(ax, 'Time [s]', 'Interpreter','latex','FontName','Times New Roman');
        end

        % styling
        grid(ax,'on'); grid(ax,'minor');
        set(ax,'FontName','Times New Roman','TickLabelInterpreter','latex', 'FontSize', 22);
    end
end

% save outputs
out_base = fullfile(results_dir, sprintf('%s_phase_errors_zoom','csm'));
exportgraphics(fig, [out_base,'.pdf'], 'ContentType','vector');
savefig(fig, [out_base,'.fig']);

%% Auxiliary functions
function [rx_signal_model_inputs, gen_kf_cfg, init_estimates_csm, init_estimates_cpssm, init_estimates_none, ar_phase_idx] = ...
    get_overall_cfgs(cache_dir, scint_model, severity, is_refractive_effects_removed_received_signal, sigma2_W_3, sampling_interval, settling_time, simulation_time, seed)
    % Define overall settings for the KF framework setup

    % General parameters
    doppler_profile = [0, 1000, 0.94];
    L1_C_over_N0_dBHz = 42;

    % Parameters for training the AR models for scintillation phase
    training_simulation_time = 300;
    is_refractive_effects_removed_training_data = true; % Exclude the refractive effects
    is_unwrapping_used = false; % This flag forces to use the wrapped phase for training the AR model

    % Parts for building the received signal
    csm_first_part = {L1_C_over_N0_dBHz, 'CSM', doppler_profile};
    csm_second_part = {'simulation_time', simulation_time, 'settling_time', settling_time};
    cpssm_first_part = {L1_C_over_N0_dBHz, 'TPPSM', doppler_profile};
    cpssm_second_part = {'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', is_refractive_effects_removed_received_signal};

    % CSM parameters for weak and strong scintillation
    S4_preset = [0.2, 0.9];
    tau0_preset = [1, 0.2];

    % CPSSM timing scaling parameters for weak and strong scintillation
    rhof_veff_ratio_preset = [1.5, 0.27]; % See right plot of Table 3.2 of my dissertation.

    switch severity
        case 'weak'
            train_cfg_csm = struct('scintillation_model', 'CSM', 'S4', S4_preset(1), 'tau0', tau0_preset(1), ...
                'simulation_time', training_simulation_time, ...
                'sampling_interval', sampling_interval, ...
                'is_unwrapping_used', is_unwrapping_used);
            train_cfg_cpssm  = struct('scintillation_model', 'TPPSM', 'scenario', severity, ...
                'rhof_veff_ratio', rhof_veff_ratio_preset(1), ...
                'simulation_time', training_simulation_time, ...
                'is_refractive_effects_removed', is_refractive_effects_removed_training_data, ...
                'sampling_interval', sampling_interval, ...
                'is_unwrapping_used', is_unwrapping_used, ...
                'is_enable_cmd_print', false);
            switch scint_model
                case 'csm'
                    ar_model_order = 6; % See right plot of Figure 4.2 of my dissertation.
                    rx_signal_model_inputs = [csm_first_part(:)',{seed},{'S4'},{S4_preset(1)},{'tau0'},{tau0_preset(1)},csm_second_part(:)'];
                case 'cpssm'
                    ar_model_order = 14; % See right plot of Figure 4.7 of my dissertation.
                    rx_signal_model_inputs = [cpssm_first_part(:)',{seed},{'tppsm_scenario'}, {'weak'},cpssm_second_part(:)', 'rhof_veff_ratio', rhof_veff_ratio_preset(1)];
            end
        case 'strong'
            train_cfg_csm = struct('scintillation_model', 'CSM', 'S4', S4_preset(2), 'tau0', tau0_preset(2), ...
                'simulation_time', training_simulation_time, ...
                'sampling_interval', sampling_interval, ...
                'is_unwrapping_used', is_unwrapping_used);
            train_cfg_cpssm  = struct('scintillation_model', 'TPPSM', 'scenario', severity, ...
                'rhof_veff_ratio', rhof_veff_ratio_preset(2), ...
                'simulation_time', training_simulation_time, ...
                'is_refractive_effects_removed', is_refractive_effects_removed_training_data, ...
                'sampling_interval', sampling_interval, ...
                'is_unwrapping_used', is_unwrapping_used, ...
                'is_enable_cmd_print', false);
            switch scint_model
                case 'csm'
                    rx_signal_model_inputs = [csm_first_part(:)',{seed},{'S4'},{S4_preset(2)},{'tau0'},{tau0_preset(2)},csm_second_part(:)'];
                    ar_model_order = 5; % See right plot of Figure 4.2 of my dissertation.
                case 'cpssm' 
                    ar_model_order = 1; % See right plot of Figure 4.7 of my dissertation.
                    rx_signal_model_inputs = [cpssm_first_part(:)',{seed},{'tppsm_scenario'},{'strong'},cpssm_second_part(:)', 'rhof_veff_ratio', rhof_veff_ratio_preset(2)];
            end
    end

    train_cfg_none  = struct('scintillation_model', 'none', 'sampling_interval', sampling_interval);

    expected_doppler_profile = [0,1000,0.94];

    gen_cfg_csm = struct( ...
      'kf_type', 'standard', ...
      'discrete_wiener_model_config', { {1, 3, sampling_interval, [0, 0, sigma2_W_3], 1} }, ...
      'scintillation_training_data_config', train_cfg_csm, ...
      'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
      'initial_states_distributions_boundaries', { {[-pi, pi], [-25, 25], [-1, 1]} }, ...
      'expected_doppler_profile', expected_doppler_profile, ...
      'augmentation_model_initializer', struct('id', 'aryule', 'model_params', struct('model_order', ar_model_order)), ...
      'is_use_cached_settings', false, ...
      'is_generate_random_initial_estimates', true, ...
      'is_enable_cmd_print', false ...
    );

    gen_cfg_cpssm = gen_cfg_csm;
    gen_cfg_cpssm.scintillation_training_data_config = train_cfg_cpssm;
    
    gen_cfg_none = gen_cfg_csm;
    gen_cfg_none.scintillation_training_data_config = train_cfg_none;
    gen_cfg_none.augmentation_model_initializer.id = 'none';
    gen_cfg_none.augmentation_model_initializer.model_params = struct();
    is_enable_cmd_print = false;

    [~, init_estimates_csm] = get_kalman_pll_config(gen_cfg_csm, cache_dir, is_enable_cmd_print);
    [~, init_estimates_cpssm] = get_kalman_pll_config(gen_cfg_cpssm, cache_dir, is_enable_cmd_print);
    [gen_kf_cfg, init_estimates_none] = get_kalman_pll_config(gen_cfg_none, cache_dir, is_enable_cmd_print);

    ar_phase_idx = length(expected_doppler_profile) + 1;
end

function [ahl_kf_ar_cfg, online_mdl_learning_cfg] = get_adaptive_cfgs()
    % Define approaches settings
    sampling_interval = 1e-2;

    % Hard-Limiting constraint threshold.
    lambda = 35;  
    
    % NWPR parameters
    T_bit = 1/50;
    M_nwpr = T_bit / sampling_interval;
    N_nwpr = 20;
    
    % For AR (KFAR) with hard limiting enabled (AHL-KF):
    ahl_kf_ar_cfg = struct(...
        'measurement_cov_adapt_algorithm', 'nwpr', ...
        'measurement_cov_adapt_algorithm_params', struct('N_nwpr', N_nwpr, 'M_nwpr', M_nwpr), ...
        'states_cov_adapt_algorithm', 'none', ...
        'hl_start_time', 50, ...
        'sampling_interval', sampling_interval, ...
        'hard_limited', struct('is_used', true, 'L1_C_over_N0_dBHz_threshold', lambda));
    
    % Online model learning setting
    online_mdl_learning_cfg = struct('is_online', false);
end
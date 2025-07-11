% error_series_ahl_kf_ar_vs_sigma2_W.m
%
% Plot AHL-KF-AR phase error time series for varying σ²_{W,3}
% in 3×2 tiled figures (rows=phase component, cols=severity) for CSM,
% CPSSM without refractive, and CPSSM with refractive.
% Mirrors error_series_comparison.m structure but loops σ² values for AHL-KF-AR only.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;
addpath(genpath(fullfile(pwd,'..','..','..','libs')));

%% 1) Reproducibility & Config
seed              = 6;
rng(seed);
sampling_interval = 1e-2;
settling_time     = 50;
simulation_time   = 300;
cache_dir         = fullfile(fileparts(mfilename('fullpath')),'cache');

%% 2) Zoom window
zoom_start = 50;  % seconds
zoom_end   = 150; % seconds
idx_full   = round(zoom_start/sampling_interval):round(zoom_end/sampling_interval);
time_zoom  = (zoom_start:sampling_interval:zoom_end).';

%% 3) AHL-KF-AR config
[ahl_kf_ar_cfg, online_mdl_learning_cfg] = get_adaptive_cfgs();

%% 4) AR-noise sweep
sigma2_W_3_sweep = logspace(-2, -10, 5);
colors           = cool(numel(sigma2_W_3_sweep));
labels_sigma     = arrayfun(@(s) sprintf('$\\sigma^2_{W,3}=%.1e$',s), sigma2_W_3_sweep, 'UniformOutput', false);

%% 5) Layout & labels
model_list   = {'csm', 'cpssm_wo_refr', 'cpssm_w_refr'};
severities   = {'weak', 'strong'};
field_labels = { ...
    '$\hat{\phi}_{\mathrm{T},1}[k\mid k-1] - \phi_{\mathrm{T},1}[k]$', ...
    '$\hat{\phi}_{\mathrm{W},1}[k\mid k-1] - \phi_{\mathrm{LOS},1}[k]$', ...
    '$\hat{\phi}_{\mathrm{AR},1}[k\mid k-1]$' ...
};
font_size    = 14;

%% 6) Loop through each model
for m = 1:numel(model_list)
    mdl = model_list{m};
    switch mdl
        case 'csm'
            scint_cfg      = 'csm';
            training_scint = 'CSM';
            refr_flag      = true;
        case 'cpssm_wo_refr'
            scint_cfg      = 'cpssm';
            training_scint = 'TPPSM';
            refr_flag      = true;
        case 'cpssm_w_refr'
            scint_cfg      = 'cpssm';
            training_scint = 'TPPSM';
            refr_flag      = false;
    end

    fig = figure('Units','normalized','Position',[0.05,0.05,0.9,0.9],'Color','w');
    tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

    for row = 1:3
        for col = 1:2
            ax = nexttile; hold(ax,'on');
            sev = severities{col};

            % Pre-get windowed true phases using first sigma
            [rx0, gen0, init_csm0, ~, ~, ar_idx] = get_overall_cfgs(...
                cache_dir, scint_cfg, sev, refr_flag, sigma2_W_3_sweep(1), ...
                sampling_interval, settling_time, simulation_time, seed);
            [sig0, true_los0, ~, diffr0, refr0] = get_received_signal(rx0{:});
            tl = true_los0(idx_full);
            df = diffr0(idx_full);
            rf = refr0(idx_full);

            % Loop sigma
            for k = 1:numel(sigma2_W_3_sweep)
                s2 = sigma2_W_3_sweep(k);
                [rx, gen, init_csm, ~, ~, ar_idx] = get_overall_cfgs(...
                    cache_dir, scint_cfg, sev, refr_flag, s2, ...
                    sampling_interval, settling_time, simulation_time, seed);
                [sig, ~, ~, ~, ~] = get_received_signal(rx{:});

                [est, ~] = get_kalman_pll_estimates(...
                    sig, gen, init_csm, 'standard', training_scint, ahl_kf_ar_cfg, online_mdl_learning_cfg);
                phi_W  = est(:,1);
                phi_AR = est(:,ar_idx);

                switch row
                    case 1  % total-error
                        err = phi_W(idx_full) + phi_AR(idx_full) - (tl + df + (~refr_flag).*rf);
                        plot(ax, time_zoom, err, 'LineWidth',1.5, 'Color',colors(k,:), 'DisplayName',labels_sigma{k});
                    case 2  % LOS-error
                        err = phi_W(idx_full) - tl;
                        plot(ax, time_zoom, err, 'LineWidth',1.5, 'Color',colors(k,:), 'DisplayName',labels_sigma{k});
                        % refractive overlay & legend dashed-dot
                        if ~refr_flag
                            h_refr = plot(ax, time_zoom, rf, '-.', 'LineWidth',1.5, 'Color','k', 'DisplayName','Refractive Phase  ($\phi_{\mathrm{R}, 1}[k]$)');
                            legend(ax, h_refr, 'Interpreter','latex','Location','best','FontName','Times New Roman');
                        end
                    case 3  % φ_AR itself
                        plot(ax, time_zoom, phi_AR(idx_full), 'LineWidth',1.5, 'Color',colors(k,:), 'DisplayName',labels_sigma{k});
                        % diffractive overlay & legend dotted
                        h_diff = plot(ax, time_zoom, df, ':', 'LineWidth',1.5, 'Color','k', 'DisplayName','Diffractive Phase ($\phi_{\mathrm{D}, 1}[k]$)');
                        legend(ax, h_diff, 'Interpreter','latex','Location','best','FontName','Times New Roman');
                end
            end

            % zero baseline only on row 1
            if row == 1
                plot(ax, time_zoom, zeros(size(idx_full)), 'k','LineWidth',2.5,'DisplayName','Target Error');
                if col == 1
                    legend(ax,'Location','best','Interpreter','latex','FontName','Times New Roman');
                end
            end

            hold(ax,'off');
            % title & labels
            if strcmp(sev, 'weak')
                title_name_final = 'Weak';
            elseif strcmp(sev, 'strong')
                title_name_final = 'Strong';
            end
            if row == 1
                title(ax, title_name_final, 'Interpreter','latex','FontName','Times New Roman');
            end
            if col == 1
                ylabel(ax, field_labels{row}, 'Interpreter','latex','FontName','Times New Roman');
            end
            if row == 3
                xlabel(ax,'Time [s]','Interpreter','latex','FontName','Times New Roman');
            end
            grid(ax,'on'); grid(ax,'minor');
            set(ax,'FontName','Times New Roman','TickLabelInterpreter','latex','FontSize',font_size);
        end
    end

    % save
    out = fullfile('results', sprintf('%s_ahl_kf_ar_vs_sigma2', mdl));
    exportgraphics(fig,[out,'.pdf'],'ContentType','vector');
    savefig(fig,[out,'.fig']);
end

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
      'discrete_wiener_model_config', { {1, 3, sampling_interval, [1e-2, 0, sigma2_W_3], 1} }, ...
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
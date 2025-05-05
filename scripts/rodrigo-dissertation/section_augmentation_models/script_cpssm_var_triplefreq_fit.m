% script_cpssm_var_triplefreq_fit.m
%
% Script to identify the goodness of fit of the VAR model for
% triple-frequency scintillation amplitude and phase time series
% genenrated by the compact phase-screen-based scintillation model (CPSSM)
% using the ARFIT algorithm [1].
%
% References:
% [1] Schneider, Tapio, and Arnold Neumaier. “Algorithm 808: ARfit—a Matlab
% Package for the Estimation of Parameters and Eigenmodes of Multivariate
% Autoregressive Models.” ACM Trans. Math. Softw. 27, no. 1
% (March 1, 2001): 58–65. https://doi.org/10.1145/382043.382316.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;
addpath(genpath(fullfile(pwd,'..','..','..', 'libs')));

% Setup output folders
fig_dir = 'pdf_figures_cpssm_triplefreq'; % where we’ll write vector PDFs
csv_dir = 'csv_data_cpssm_triplefreq';    % where we’ll write CSV tables

if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end
if ~exist(csv_dir,'dir')
    mkdir(csv_dir);
end

%% Simulation parameters

simulation_time = 300;
sampling_interval = 0.01;
severities = {'Weak','Moderate','Strong'};
frequency_bands = {'L1', 'L2', 'L5'};
freq_amount = numel(frequency_bands);
cpssm_params = struct( ...
    'Weak',    {'weak',     'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 1.5},...
    'Moderate',{'moderate', 'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 0.8},...
    'Strong',  {'strong',   'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 0.27}...
    );

%% Monte Carlo optimal VAR model order assessment
% Reduce the amount of `mc_runs` for faster results.
mc_runs    = 100;
min_order  = 1;
max_order  = 30;
optimal_orders_amp_total_phs   = zeros(mc_runs,numel(severities));
sbc_amp_total_phs_array = zeros(mc_runs, numel(severities), max_order - min_order + 1);

seed = 1;
for mc_idx = 1:mc_runs
    for i = 1:numel(severities)
        severity = severities{i};
        rng(seed);
        [scint_ts, ~] = get_tppsm_data(cpssm_params.(severity), 'seed', seed);
        amp_ts    = abs(scint_ts);

        total_phs_ts_training  = zeros(size(amp_ts));
        for freq_idx = 1:freq_amount
            % NOTE: `get_corrected_phase` is a function from the
            % 'gnss_scintillation_simulator' git submodule.
            total_phs_ts_training(:,freq_idx) = get_corrected_phase(scint_ts(:,freq_idx));
        end
        amp_total_phs_ts_training = [amp_ts, total_phs_ts_training];

        % Fit the scintillation amplitude, total, refractive and
        % diffractive phases, respectively to a AR model, respectively.
        [~, A_amp_total_phs, ~, sbc_amp_total_phs]   = arfit(amp_total_phs_ts_training,   min_order, max_order);

        optimal_orders_amp_total_phs(mc_idx,i)   = size(A_amp_total_phs,2)/size(A_amp_total_phs,1);

        sbc_amp_total_phs_array(mc_idx, i, :) = sbc_amp_total_phs;

        seed = seed + 1;
    end
end

orders       = min_order:max_order;
counts_amp_phs   = zeros(numel(orders),numel(severities));
for i = 1:numel(severities)
    counts_amp_phs(:,i)       = histcounts(optimal_orders_amp_total_phs(:,i),   [orders, orders(end)+1]);
end

% Normalize counts to percentages
pct_amp_total_phs = counts_amp_phs / mc_runs * 100;

%% Plot model order selection frequency

%figure('Position',[50,50,1100,400]);
figure;
colors = lines(numel(severities));  % or get(gca,'ColorOrder')

% Amplitude
hold on;
for j = 1:numel(severities)
    plot(orders, pct_amp_total_phs(:,j), '-o', ...
        'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
        'DisplayName',severities{j});
end
hold off;
xlabel('VAR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal VAR Order – Amplitudes and total phases');
legend('Location','best');
grid on;

%% Export Optimal VAR Order Plot and CSV

% 1) Export the current figure as a cropped, fully-vector PDF
fig_name = 'optimal_var_order_pct';
pdf_file = fullfile(fig_dir, [fig_name, '.pdf']);
exportgraphics(gcf, pdf_file, 'ContentType', 'vector');

% 2) Build and export CSV of percentage data for later TikZ plotting
% Column order: Order, Amp_Weak, Amp_Moderate, Amp_Strong,
%               Tot_Weak, Tot_Moderate, Tot_Strong,
%               Refr_Weak, Refr_Moderate, Refr_Strong,
%               Diff_Weak, Diff_Moderate, Diff_Strong

T = table( orders.', ...
    pct_amp_total_phs(:,1), pct_amp_total_phs(:,2), pct_amp_total_phs(:,3), ...
    'VariableNames', { ...
    'Order','Amp_total_phs_Weak', ...
    'Amp_total_phs_Moderate','Amp_total_phs_Strong'} );

csv_file = fullfile(csv_dir, [fig_name, '.csv']);
writetable(T, csv_file);

%% Plot Mean SBC for Amplitudes & Phase Components
orders               = min_order:max_order;
mean_sbc_amp_total_phs         = squeeze(mean(sbc_amp_total_phs_array, 1)).';   % [severity × order]

%figure('Position',[50,50,1100,400]);
figure;
colors = get(gca,'ColorOrder');
plot(orders, mean_sbc_amp_total_phs, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_amp_total_phs,2)
    [minVal, minIdx] = min(mean_sbc_amp_total_phs(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('VAR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Amplitude');
legend(severities,'Location','best'); grid on; hold off;

%% Export Mean SBC Plot and CSV

% Export Mean SBC figure as a cropped, fully-vector PDF
fig_name = 'mean_sbc';
pdf_file = fullfile(fig_dir, [fig_name, '.pdf']);
exportgraphics(gcf, pdf_file, 'ContentType','vector');

% Build and export CSV of Mean SBC data for TikZ
T_sbc = table(orders.', ...
    mean_sbc_amp_total_phs(:,1),    mean_sbc_amp_total_phs(:,2),    mean_sbc_amp_total_phs(:,3), ...
    'VariableNames',{ ...
    'Order', ...
    'Amp_total_phs_Weak','Amp_total_phs_Moderate','Amp_total_phs_Strong'} );
csv_file = fullfile(csv_dir, [fig_name,'.csv']);
writetable(T_sbc, csv_file);

%% Residual Analysis
[~, highest_freq_idx_amp]   = max(counts_amp_phs,[],1);
[~, highest_freq_idx_total_phs] = max(counts_total_phs,[],1);
[~, highest_freq_idx_refr_phs] = max(counts_refr_phs,[],1);
[~, highest_freq_idx_diff_phs] = max(counts_diff_phs,[],1);

residuals = struct('amplitude',[] ,'total_phs', [], 'refr_phs', [], 'diff_phs', []);
residuals_downsamp = struct('amplitude',[] ,'total_phs', [], 'refr_phs', [], 'diff_phs', []);

seed = 10;
for i = 1:numel(severities)
    severity = severities{i};
    rng(seed);

    % Training the VAR model
    [scint_ts_training, refr_phs_ts_training] = get_tppsm_data(cpssm_params.(severity), 'seed', seed);
    amp_ts_training = abs(scint_ts);

    total_phs_ts_training = zeros(size(amp_ts_training));
    for freq_idx = 1:freq_amount
        total_phs_ts_training(:,freq_idx) = ...
            get_corrected_phase(scint_ts_training(:,freq_idx));
    end

    diff_phs_ts_training = wrapToPi(total_phs_ts_training - refr_phs_ts_training);

    ord_amp = orders(highest_freq_idx_amp(i));
    ord_total_phs = orders(highest_freq_idx_total_phs(i));
    ord_refr_phs = orders(highest_freq_idx_refr_phs(i));
    ord_diff_phase = orders(highest_freq_idx_diff_phs(i));

    [w_amp, A_amp_total_phs]   = arfit(amp_ts_training, ord_amp, ord_amp);
    [w_total_phs, A_total_phs] = arfit(total_phs_ts_training, ord_total_phs, ord_total_phs);
    [w_refr_phs, A_refr_phs] = arfit(refr_phs_ts_training, ord_refr_phs, ord_refr_phs);
    [w_diff_phase, A_diff_phase] = arfit(diff_phs_ts_training, ord_diff_phase, ord_diff_phase);

    % Obtaining the residuals
    % NOTE: Another generation seed is used herein to evaluate the reisuduals
    % correlation with statistical significance.
    rng(seed+1);
    [scint_ts_for_residue, refr_phs_ts_training] = get_tppsm_data(cpssm_params.(severity), 'seed', seed + 1);
    amp_ts_for_residue = abs(scint_ts_for_residue);

    total_phs_ts_for_residue = zeros(size(amp_ts_for_residue));
    for freq_idx = 1:freq_amount
        total_phs_ts_for_residue(:,freq_idx) = ...
            get_corrected_phase(scint_ts_for_residue(:,freq_idx));
    end
    diff_phs_ts_for_residue = wrapToPi(total_phs_ts_for_residue - refr_phs_ts_training);

    [~, res_amp]  = arres(w_amp, A_amp_total_phs, amp_ts_for_residue, 60);
    [~, res_total_phs] = arres(w_total_phs, A_total_phs, total_phs_ts_for_residue, 60);
    [~, res_refr_phs] = arres(w_refr_phs, A_refr_phs, refr_phs_ts_training, 60);
    [~, res_diff_phase] = arres(w_diff_phase, A_diff_phase, diff_phs_ts_for_residue, 60);

    residuals.amplitude.(severity) = [NaN(ord_amp, size(res_amp,2)); res_amp];
    residuals.total_phs.(severity) = [NaN(ord_total_phs, size(res_total_phs,2));   res_total_phs];
    residuals.refr_phs.(severity) = [NaN(ord_refr_phs, size(res_refr_phs,2));   res_refr_phs];
    residuals.diff_phs.(severity) = [NaN(ord_diff_phase, size(res_diff_phase,2));   res_diff_phase];

    seed = seed + 1;
end

%% Plot the downsampled residuals of the fitted model

time = sampling_interval:sampling_interval:simulation_time;
% Downsample the time series for lower file size, if necessary.
downsamp_factor = 1;
downsamp_sampling_interval = downsamp_factor*sampling_interval;
time_downsamp = downsample(time, downsamp_factor);

for i = 1:numel(severities)
    severity = severities{i};
    residuals_downsamp.amplitude.(severity) = downsample(residuals.amplitude.(severity),downsamp_factor);
    residuals_downsamp.total_phs.(severity) = downsample(residuals.total_phs.(severity),downsamp_factor);
    residuals_downsamp.refr_phs.(severity) = downsample(residuals.refr_phs.(severity),downsamp_factor);
    residuals_downsamp.diff_phs.(severity) = downsample(residuals.diff_phs.(severity),downsamp_factor);
end

% Define plot order and time axis
plot_order = {'Strong','Moderate','Weak'};
% Invert default colormap so
% strong→yellow (3rd), moderate→red (2nd), weak→blue (1st)
base_colors = lines(numel(plot_order));
colors      = base_colors([3,2,1],:);

for freq_idx = 1:numel(frequency_bands)
    figure('Position',[50,50,1100,600]);

    % Amplitude Residuals
    subplot(2,2,1); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        selected_amp_res = residuals_downsamp.amplitude.(sev);
        plot(time_downsamp, selected_amp_res(:,freq_idx), ...
            'LineWidth',1, 'Color',colors(k,:), ...
            'DisplayName', sev);
    end
    hold off;
    xlabel('Time [s]'); ylabel('Residuals');
    title(['Amplitude Residuals - ', frequency_bands{freq_idx}]);
    legend({'Strong','Moderate','Weak'},'Location','best');
    grid on;

    % Total Phase Residuals
    subplot(2,2,2); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        selected_total_phs = residuals_downsamp.total_phs.(sev);
        plot(time_downsamp, selected_total_phs(:,freq_idx), ...
            'LineWidth',1, 'Color',colors(k,:), ...
            'DisplayName', sev);
    end
    hold off;
    xlabel('Time [s]'); ylabel('Residuals');
    title(['Total Phase Residuals - ', frequency_bands{freq_idx}]);
    legend({'Strong','Moderate','Weak'},'Location','best');
    grid on;

    % Refractive Phase Residuals
    subplot(2,2,3); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        selected_refr_phs = residuals_downsamp.refr_phs.(sev);
        plot(time_downsamp, selected_refr_phs(:,freq_idx), ...
            'LineWidth',1, 'Color',colors(k,:), ...
            'DisplayName', sev);
    end
    hold off;
    xlabel('Time [s]'); ylabel('Residuals');
    title(['Refractive Phase Residuals - ', frequency_bands{freq_idx}]);
    legend({'Strong','Moderate','Weak'},'Location','best');
    grid on;

    % 4) Diffractive Phase Residuals
    subplot(2,2,4); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        selected_diff_phs = residuals_downsamp.diff_phs.(sev);
        plot(time_downsamp, selected_diff_phs(:,freq_idx), ...
            'LineWidth',1, 'Color',colors(k,:), ...
            'DisplayName', sev);
    end
    hold off;
    xlabel('Time [s]'); ylabel('Residuals');
    title(['Diffractive Phase Residuals - ', frequency_bands{freq_idx}]);
    legend({'Strong','Moderate','Weak'},'Location','best');
    grid on;

    %% Export Residuals Plot and CSV
    % Export residuals figure as a cropped, fully-vector PDF
    fig_name = ['residuals', frequency_bands{freq_idx}];
    pdf_file = fullfile(fig_dir, [fig_name, '.pdf']);
    exportgraphics(gcf, pdf_file, 'ContentType','vector');

    % Build and export CSV of residuals data for TikZ
    T = table(time_downsamp.', ...
        residuals.amplitude.Weak,   residuals.amplitude.Moderate,   residuals.amplitude.Strong, ...
        residuals.total_phs.Weak,   residuals.total_phs.Moderate,   residuals.total_phs.Strong, ...
        residuals.refr_phs.Weak,    residuals.refr_phs.Moderate,    residuals.refr_phs.Strong, ...
        residuals.diff_phs.Weak,  residuals.diff_phs.Moderate,  residuals.diff_phs.Strong, ...
        'VariableNames',{ ...
        'Time_s', ...
        'Amp_Weak','Amp_Moderate','Amp_Strong', ...
        'Tot_Weak','Tot_Moderate','Tot_Strong', ...
        'Refr_Weak','Refr_Moderate','Refr_Strong', ...
        'Diff_Weak','Diff_Moderate','Diff_Strong' } );
    csv_file = fullfile(csv_dir, [fig_name, '.csv']);
    writetable(T, csv_file);

    %% Compute one-sided residuals ACFs

    % Amount of lags on the ACF
    lags = 20;
    acfs = struct('amplitude',[] ,'total_phs', [], 'refr_phs', [], 'diff_phs', []);

    for i = 1:numel(severities)
        sev = severities{i};

        % Removing the NaNs of the amplitude and phase residuals
        amp_res  = residuals.amplitude.(severity);
        amp_res  = amp_res(~isnan(amp_res));
        total_phs_res  = residuals.total_phs.(severity);
        total_phs_res  = total_phs_res(~isnan(total_phs_res));
        refr_phs_res  = residuals.refr_phs.(severity);
        refr_phs_res  = refr_phs_res(~isnan(refr_phs_res));
        diff_phase_res  = residuals.diff_phs.(severity);
        diff_phase_res  = diff_phase_res(~isnan(diff_phase_res));

        acf_amp = xcorr(amp_res, amp_res, lags, 'normalized');
        acf_total_phs = xcorr(total_phs_res, total_phs_res, lags, 'normalized');
        acf_refr_phs = xcorr(refr_phs_res, refr_phs_res, lags, 'normalized');
        acf_diff_phase = xcorr(diff_phase_res, diff_phase_res, lags, 'normalized');

        acfs.amplitude.(severity) = acf_amp(lags+1:end);
        acfs.total_phs.(severity) = acf_total_phs(lags+1:end);
        acfs.refr_phs.(severity) = acf_refr_phs(lags+1:end);
        acfs.diff_phs.(severity) = acf_diff_phase(lags+1:end);
    end
 
    plot_order = {'Strong','Moderate','Weak'};
    time_lag   = (0:lags) * sampling_interval;
    colors     = lines(numel(plot_order));
    colors     = colors([3,2,1],:);    % strong→yellow, moderate→red, weak→blue
    stem_width = 1.5;
    markers    = struct('Weak','o','Moderate','s','Strong','^');

    % Plot ACFs in 2×2 grid: amplitude, total phase, refractive, diffractive
    figure('Position',[50,50,1100,600]);

    % Amplitude ACF
    subplot(2,2,1); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        stem(time_lag, acfs.amplitude.(sev), 'LineWidth',stem_width, ...
            'Color',colors(k,:), 'Marker',markers.(sev), ...
            'DisplayName',sev);
    end
    hold off;
    xlabel('Time Lag [s]'); ylabel('Normalized ACF');
    title('Amplitude Residuals ACF');
    legend('Location','best'); grid on;

    % Total Phase ACF
    subplot(2,2,2); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        stem(time_lag, acfs.total_phs.(sev), 'LineWidth',stem_width, ...
            'Color',colors(k,:), 'Marker',markers.(sev), ...
            'DisplayName',sev);
    end
    hold off;
    xlabel('Time Lag [s]'); ylabel('Normalized ACF');
    title('Total Phase Residuals ACF');
    legend('Location','best'); grid on;

    % Refractive Phase ACF
    subplot(2,2,3); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        stem(time_lag, acfs.refr_phs.(sev), 'LineWidth',stem_width, ...
            'Color',colors(k,:), 'Marker',markers.(sev), ...
            'DisplayName',sev);
    end
    hold off;
    xlabel('Time Lag [s]'); ylabel('Normalized ACF');
    title('Refractive Phase Residuals ACF');
    legend('Location','best'); grid on;

    % Diffractive Phase ACF
    subplot(2,2,4); hold on;
    for k = 1:numel(plot_order)
        sev = plot_order{k};
        stem(time_lag, acfs.diff_phs.(sev), 'LineWidth',stem_width, ...
            'Color',colors(k,:), 'Marker',markers.(sev), ...
            'DisplayName',sev);
    end
    hold off;
    xlabel('Time Lag [s]'); ylabel('Normalized ACF');
    title('Diffractive Phase Residuals ACF');
    legend('Location','best'); grid on;

    %% Export the residuals ACFs plots and CSV
    % Export figure as cropped, fully‐vector PDF
    fig_name = 'residuals_acf';  % ← change this to something meaningful
    pdf_file = fullfile(fig_dir, [fig_name,'.pdf']);

    % exportgraphics crops to the tight content box and preserves vectors
    exportgraphics(gcf, pdf_file, 'ContentType','vector');

    % Export the data behind the plot as CSV

    T = table( time_lag.', ...
        acfs.amplitude.Strong,   acfs.amplitude.Moderate,   acfs.amplitude.Weak, ...
        acfs.total_phs.Strong,   acfs.total_phs.Moderate,   acfs.total_phs.Weak, ...
        acfs.refr_phs.Strong,    acfs.refr_phs.Moderate,    acfs.refr_phs.Weak, ...
        acfs.diff_phs.Strong,  acfs.diff_phs.Moderate,  acfs.diff_phs.Weak, ...
        'VariableNames',{ ...
        'Lag_s', ...
        'Amp_Strong','Amp_Moderate','Amp_Weak', ...
        'Tot_Strong','Tot_Moderate','Tot_Weak', ...
        'Refr_Strong','Refr_Moderate','Refr_Weak', ...
        'Diff_Strong','Diff_Moderate','Diff_Weak' } );

    writetable(T, fullfile(csv_dir, [fig_name,'.csv']));
end
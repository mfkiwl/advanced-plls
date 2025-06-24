% script_cpssm_ar_singlefreq_fit.m
%
% Script to identify the goodness of fit of the AR model for
% single-frequency scintillation amplitude and phase time series
% genenrated by the compact phase-screen-based scintillation model (CPSSM)
% using the ARFIT algorithm [1].
%
% References:
% [1] Schneider, Tapio, and Arnold Neumaier. “Algorithm 808: ARfit—a Matlab
% Package for the Estimation of Parameters and Eigenmodes of Multiariate
% Autoregressive Models.” ACM Trans. Math. Softw. 27, no. 1
% (March 1, 2001): 58–65. https://doi.org/10.1145/382043.382316.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;
addpath(genpath(fullfile(pwd,'..','..','..', 'libs')));

% Setup output folders
fig_dir = 'pdf_figures_cpssm_singlefreq'; % where we’ll write vector PDFs
csv_dir = 'csv_data_cpssm_singlefreq';    % where we’ll write CSV tables

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
cpssm_params = struct( ...
    'Weak',    {'weak',     'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 1.5},...
    'Moderate',{'moderate', 'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 0.8},...
    'Strong',  {'strong',   'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 0.27}...
    );
font_size = 11;

%% Monte Carlo optimal AR model order assessment
mc_runs    = 300;
min_order  = 1;
max_order  = 30;
optimal_orders_amp   = zeros(mc_runs,numel(severities));
optimal_orders_total_phs = zeros(mc_runs,numel(severities));
optimal_orders_refr_phs = zeros(mc_runs,numel(severities));
optimal_orders_diff_phs = zeros(mc_runs,numel(severities));
sbc_amp_array = zeros(mc_runs, numel(severities), max_order - min_order + 1);
sbc_total_phs_array = zeros(mc_runs, numel(severities), max_order - min_order + 1);
sbc_refr_phs_array = zeros(mc_runs, numel(severities), max_order - min_order + 1);
sbc_diff_phs_array = zeros(mc_runs, numel(severities), max_order - min_order + 1);

seed = 1;
for mc_idx = 1:mc_runs
    for i = 1:numel(severities)
        severity = severities{i};
        rng(seed);
        [scint_ts, refr_phs_ts_training] = get_tppsm_multifreq_data(cpssm_params.(severity), 'seed', seed);
        amp_ts    = abs(scint_ts(:,1));
        % NOTE: `get_corrected_phase` is a function from the
        % 'gnss_scintillation_simulator' git submodule.
        total_phs_ts_training  = get_corrected_phase(scint_ts(:,1));
        % Get the diffractive phase as the wrapped version of the
        % difference between the propagated field's phase and the
        % refractive phase time series at the ionospheric piercing point
        % (IPP).
        % NOTE: The `wrapToPi` function serves to wrap the diffractive phase
        % time series within  the [-pi, pi] bounds, in order to be
        % comparable with the CSM model in a feasible way.
        diff_phs_ts_training = wrapToPi(total_phs_ts_training - refr_phs_ts_training(:,1));

        % Fit the scintillation amplitude, total, refractive and
        % diffractive phases, respectively to a AR model, respectively.
        [~, A_amp, ~, sbc_amp]   = arfit(amp_ts,   min_order, max_order);
        [~, A_total_phs, ~, sbc_total_phs] = arfit(total_phs_ts_training, min_order, max_order);
        [~, A_refr_phs, ~, sbc_refr_phs] = arfit(refr_phs_ts_training(:,1), min_order, max_order);
        [~, A_diff_phs, ~, sbc_diff_phs] = arfit(diff_phs_ts_training, min_order, max_order);

        optimal_orders_amp(mc_idx,i)   = size(A_amp,2)/size(A_amp,1);
        optimal_orders_total_phs(mc_idx,i) = size(A_total_phs,2)/size(A_total_phs,1);
        optimal_orders_refr_phs(mc_idx,i) = size(A_refr_phs,2)/size(A_refr_phs,1);
        optimal_orders_diff_phs(mc_idx,i) = size(A_diff_phs,2)/size(A_diff_phs,1);

        sbc_amp_array(mc_idx, i, :) = sbc_amp;
        sbc_total_phs_array(mc_idx, i, :) = sbc_total_phs;
        sbc_refr_phs_array(mc_idx, i, :) = sbc_refr_phs;
        sbc_diff_phs_array(mc_idx, i, :) = sbc_diff_phs;

        seed = seed + 1;
    end
end

orders       = min_order:max_order;
counts_amp   = zeros(numel(orders),numel(severities));
counts_total_phs = zeros(numel(orders),numel(severities));
counts_refr_phs = zeros(numel(orders),numel(severities));
counts_diff_phs = zeros(numel(orders),numel(severities));
for i = 1:numel(severities)
    counts_amp(:,i)   = histcounts(optimal_orders_amp(:,i),   [orders, orders(end)+1]);
    counts_total_phs(:,i) = histcounts(optimal_orders_total_phs(:,i), [orders, orders(end)+1]);
    counts_refr_phs(:,i) = histcounts(optimal_orders_refr_phs(:,i), [orders, orders(end)+1]);
    counts_diff_phs(:,i) = histcounts(optimal_orders_diff_phs(:,i), [orders, orders(end)+1]);
end

% Normalize counts to percentages
pct_amp       = counts_amp       / mc_runs * 100;
pct_total_phs = counts_total_phs / mc_runs * 100;
pct_refr_phs  = counts_refr_phs  / mc_runs * 100;
pct_diff_phs  = counts_diff_phs  / mc_runs * 100;

%% Plot model order selection frequency
figure('Position',[100,100,1000,500]);
colors = lines(numel(severities));  % or get(gca,'ColorOrder')
markers = {'o', 's', '^'};

% Amplitude
subplot(2,2,1);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_amp(:,j), ['-', markers{j}], ...
        'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
        'DisplayName',severities{j});
end
hold off;
xlabel('AR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal AR Order – Amplitude');
legend('Location','best');
set(gca, 'FontSize', font_size);
grid on;

% Total phase
subplot(2,2,2);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_total_phs(:,j), ['-', markers{j}], ...
        'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
        'DisplayName',severities{j});
end
hold off;
xlabel('AR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal AR Order – Total Phase');
% legend('Location','best');
set(gca, 'FontSize', font_size);
grid on;

% Refractive phase
subplot(2,2,3);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_refr_phs(:,j), ['-', markers{j}], ...
        'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
        'DisplayName',severities{j});
end
hold off;
xlabel('AR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal AR Order – Refractive Phase');
% legend('Location','best');
set(gca, 'FontSize', font_size);
grid on;

% Diffractive phase
subplot(2,2,4);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_diff_phs(:,j), ['-', markers{j}], ...
        'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
        'DisplayName',severities{j});
end
hold off;
xlabel('AR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal AR Order – Diffractive Phase');
% legend('Location','best');
set(gca, 'FontSize', font_size);
grid on;

%% Export Optimal AR Order Plot and CSV

% Export the current figure as a cropped, fully-vector PDF
fig_name = 'optimal_ar_order_frequency_singlefreq_cpssm';
pdf_file = fullfile(fig_dir, [fig_name, '.pdf']);
exportgraphics(gcf, pdf_file, 'ContentType', 'vector');

% Build and export CSV of percentage data for later TikZ plotting
% Column order: Order, Amp_Weak, Amp_Moderate, Amp_Strong,
%               Tot_Weak, Tot_Moderate, Tot_Strong,
%               Refr_Weak, Refr_Moderate, Refr_Strong,
%               Diff_Weak, Diff_Moderate, Diff_Strong

T = table( orders.', ...
    pct_amp(:,1),       pct_amp(:,2),       pct_amp(:,3), ...
    pct_total_phs(:,1), pct_total_phs(:,2), pct_total_phs(:,3), ...
    pct_refr_phs(:,1),  pct_refr_phs(:,2),  pct_refr_phs(:,3), ...
    pct_diff_phs(:,1),  pct_diff_phs(:,2),  pct_diff_phs(:,3), ...
    'VariableNames', { ...
    'Order', ...
    'Amp_Weak','Amp_Moderate','Amp_Strong', ...
    'Tot_Weak','Tot_Moderate','Tot_Strong', ...
    'Refr_Weak','Refr_Moderate','Refr_Strong', ...
    'Diff_Weak','Diff_Moderate','Diff_Strong' } );

csv_file = fullfile(csv_dir, [fig_name, '.csv']);
writetable(T, csv_file);

%% Plot Mean SBC for Amplitude & Phase Components
orders               = min_order:max_order;
mean_sbc_amp         = squeeze(mean(sbc_amp_array,         1)).';
mean_sbc_total_phs   = squeeze(mean(sbc_total_phs_array,   1)).';
mean_sbc_refr_phs    = squeeze(mean(sbc_refr_phs_array,    1)).';
mean_sbc_diff_phs    = squeeze(mean(sbc_diff_phs_array,    1)).';

figure('Position',[100,100,1000,500]);

% Amplitude
subplot(2,2,1);
colors = get(gca,'ColorOrder');
plot(orders, mean_sbc_amp, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_amp,2)
    [minVal, minIdx] = min(mean_sbc_amp(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('AR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Amplitude');
legend(severities,'Location','northeast'); grid on; hold off;
set(gca, 'FontSize', font_size);

% Total phase
subplot(2,2,2);
plot(orders, mean_sbc_total_phs, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_total_phs,2)
    [minVal, minIdx] = min(mean_sbc_total_phs(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('AR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Total Phase');
% legend(severities,'Location','best'); 
grid on; hold off;
set(gca, 'FontSize', font_size);

% Refractive phase
subplot(2,2,3);
plot(orders, mean_sbc_refr_phs, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_refr_phs,2)
    [minVal, minIdx] = min(mean_sbc_refr_phs(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('AR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Refractive Phase');
% legend(severities,'Location','best'); 
grid on; hold off;
set(gca, 'FontSize', font_size);

% Diffractive phase
subplot(2,2,4);
plot(orders, mean_sbc_diff_phs, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_diff_phs,2)
    [minVal, minIdx] = min(mean_sbc_diff_phs(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('AR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Diffractive Phase');
% legend(severities,'Location','best'); 
grid on; hold off;
set(gca, 'FontSize', font_size);

%% Export Mean SBC Plot and CSV
% Export Mean SBC figure as a cropped, fully-vector PDF
fig_name = 'mean_sbc_singlefreq_cpssm';
pdf_file = fullfile(fig_dir, [fig_name, '.pdf']);
exportgraphics(gcf, pdf_file, 'ContentType','vector');

% Build and export CSV of Mean SBC data for TikZ
T_sbc = table(orders.', ...
    mean_sbc_amp(:,1),    mean_sbc_amp(:,2),    mean_sbc_amp(:,3), ...
    mean_sbc_total_phs(:,1), mean_sbc_total_phs(:,2), mean_sbc_total_phs(:,3), ...
    mean_sbc_refr_phs(:,1),  mean_sbc_refr_phs(:,2),  mean_sbc_refr_phs(:,3), ...
    mean_sbc_diff_phs(:,1),  mean_sbc_diff_phs(:,2),  mean_sbc_diff_phs(:,3), ...
    'VariableNames',{ ...
    'Order', ...
    'Amp_Weak','Amp_Moderate','Amp_Strong', ...
    'Tot_Weak','Tot_Moderate','Tot_Strong', ...
    'Refr_Weak','Refr_Moderate','Refr_Strong', ...
    'Diff_Weak','Diff_Moderate','Diff_Strong' } );
csv_file = fullfile(csv_dir, ['mean_sbc','.csv']);
writetable(T_sbc, csv_file);

%% Residual Analysis
[~, min_sbc_idx_amp]   = min(mean_sbc_amp,[],1);
[~, min_sbc_idx_total_phs] = min(mean_sbc_total_phs,[],1);
[~, min_sbc_idx_refr_phs] = min(mean_sbc_refr_phs,[],1);
[~, min_sbc_idx_diff_phs] = min(mean_sbc_diff_phs,[],1);

residuals = struct('amplitude',[] ,'total_phs', [], 'refr_phs', [], 'diff_phs', []);
residuals_downsamp = struct('amplitude',[] ,'total_phs', [], 'refr_phs', [], 'diff_phs', []);

seed = 10;
for i = 1:numel(severities)
    severity = severities{i};
    rng(seed);

    % Training the AR model
    [scint_ts, refr_phs_ts_training] = get_tppsm_multifreq_data(cpssm_params.(severity), 'seed', seed);
    amp_ts_training    = abs(scint_ts(:,1));
    total_phs_ts_training  = get_corrected_phase(scint_ts(:,1));
    diff_phs_ts_training = wrapToPi(total_phs_ts_training - refr_phs_ts_training(:,1));

    ord_amp = orders(min_sbc_idx_amp(i));
    ord_total_phs = orders(min_sbc_idx_total_phs(i));
    ord_refr_phs = orders(min_sbc_idx_refr_phs(i));
    ord_diff_phase = orders(min_sbc_idx_diff_phs(i));

    [w_amp, A_amp]   = arfit(amp_ts_training, ord_amp, ord_amp);
    [w_total_phs, A_total_phs] = arfit(total_phs_ts_training, ord_total_phs, ord_total_phs);
    [w_refr_phs, A_refr_phs] = arfit(refr_phs_ts_training(:,1), ord_refr_phs, ord_refr_phs);
    [w_diff_phase, A_diff_phase] = arfit(diff_phs_ts_training, ord_diff_phase, ord_diff_phase);

    % Obtaining the residuals
    % NOTE: Another generation seed is used herein to evaluate the reisuduals
    % correlation with statistical significance.
    rng(seed+1);
    [scint_ts, refr_phs_ts_training] = get_tppsm_multifreq_data(cpssm_params.(severity), 'seed', seed + 1);
    amp_ts_for_residuals    = abs(scint_ts(:,1));
    total_phs_ts_for_residuals  = get_corrected_phase(scint_ts(:,1));
    diff_phs_ts_for_residuals = wrapToPi(total_phs_ts_for_residuals - refr_phs_ts_training(:,1));

    [~, res_amp]   = arres(w_amp,   A_amp,   amp_ts_for_residuals, 60);
    [~, res_total_phs]   = arres(w_total_phs,   A_total_phs,   total_phs_ts_for_residuals, 60);
    [~, res_refr_phs]   = arres(w_refr_phs,   A_refr_phs,   refr_phs_ts_training(:,1), 60);
    [~, res_diff_phase]   = arres(w_diff_phase,   A_diff_phase,   diff_phs_ts_for_residuals, 60);

    residuals.amplitude.(severity) = [NaN(ord_amp,1);   res_amp];
    residuals.total_phs.(severity) = [NaN(ord_total_phs,1);   res_total_phs];
    residuals.refr_phs.(severity) = [NaN(ord_refr_phs,1);   res_refr_phs];
    residuals.diff_phs.(severity) = [NaN(ord_diff_phase,1);   res_diff_phase];

    seed = seed + 1;
end

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

%% Plot the downsampled residuals of the fitted model

% Define plot order and time axis
plot_order = {'Strong','Moderate','Weak'};

% Invert default colormap so
% strong→yellow (3rd), moderate→red (2nd), weak→blue (1st)
base_colors = lines(numel(plot_order));
colors      = base_colors([3,2,1],:);

figure('Position',[100,100,1000,500]);

% Amplitude residuals
subplot(2,2,1); hold on;
for k = 1:numel(plot_order)
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.amplitude.(sev), ...
        'LineWidth',1, 'Color',colors(k,:), ...
        'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals');
title('Amplitude Residuals');
legend({'Strong','Moderate','Weak'},'Location','best', 'Direction', 'reverse');
set(gca, 'FontSize', font_size);
grid on;


% Total phase residuals
subplot(2,2,2); hold on;
for k = 1:numel(plot_order)
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.total_phs.(sev), ...
        'LineWidth',1, 'Color',colors(k,:), ...
        'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals [rad]');
title('Total Phase Residuals');
% legend({'Strong','Moderate','Weak'},'Location','best', 'Direction', 'reverse');
set(gca, 'FontSize', font_size);
grid on;

% Refractive Phase residuals
subplot(2,2,3); hold on;
for k = 1:numel(plot_order)
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.refr_phs.(sev), ...
        'LineWidth',1, 'Color',colors(k,:), ...
        'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals [rad]');
title('Refractive Phase Residuals');
% legend({'Strong','Moderate','Weak'},'Location','best', 'Direction', 'reverse');
set(gca, 'FontSize', font_size);
grid on;

% Diffractive Phase residuals
subplot(2,2,4); hold on;
for k = 1:numel(plot_order)
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.diff_phs.(sev), ...
        'LineWidth',1, 'Color',colors(k,:), ...
        'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals [rad]');
title('Diffractive Phase Residuals');
% legend({'Strong','Moderate','Weak'},'Location','best', 'Direction', 'reverse');
set(gca, 'FontSize', font_size);
grid on;

%% Export Residuals Plot and CSV
% Export residuals figure as a cropped, fully-vector PDF
fig_name = 'residuals_singlefreq_cpssm';
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
lags        = 20;
acfs        = struct('amplitude',[] ,'total_phs', [], 'refr_phs', [], 'diff_phs', []);

for i = 1:numel(severities)
    severity = severities{i};

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

%% Plot the residuals ACFs
plot_order = {'Strong','Moderate','Weak'};
time_lag   = (0:lags) * sampling_interval;
colors     = lines(numel(plot_order));
colors     = colors([3,2,1],:);    % strong→yellow, moderate→red, weak→blue
stem_width = 1.5;
markers    = struct('Weak','o','Moderate','s','Strong','^');

% Plot ACFs in 2×2 grid: Amplitude, total phase, refractive, diffractive
figure('Position',[100,100,1000,500]);

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
legend('Location','best', 'Direction', 'reverse'); 
grid on;
set(gca, 'FontSize', font_size);

% Total Phase ACF
subplot(2,2,2); hold on;
for k = 1:numel(plot_order)
    sev = plot_order{k};
    stem(time_lag, acfs.total_phs.(sev), 'LineWidth',stem_width, ...
        'Color',colors(k,:), 'Marker',markers.(sev), ...
        'DisplayName',sev);
end
hold off;
xlabel('Time Lag [s]'); ylabel('Normalized ACF [rad^2]');
title('Total Phase Residuals ACF');
% legend('Location','best', 'Direction', 'reverse'); 
grid on;
set(gca, 'FontSize', font_size);

% Refractive Phase ACF
subplot(2,2,3); hold on;
for k = 1:numel(plot_order)
    sev = plot_order{k};
    stem(time_lag, acfs.refr_phs.(sev), 'LineWidth',stem_width, ...
        'Color',colors(k,:), 'Marker',markers.(sev), ...
        'DisplayName',sev);
end
hold off;
xlabel('Time Lag [s]'); ylabel('Normalized ACF [rad^2]');
title('Refractive Phase Residuals ACF');
% legend('Location','best', 'Direction', 'reverse'); 
grid on;
set(gca, 'FontSize', font_size);

% 4) Diffractive Phase ACF
subplot(2,2,4); hold on;
for k = 1:numel(plot_order)
    sev = plot_order{k};
    stem(time_lag, acfs.diff_phs.(sev), 'LineWidth',stem_width, ...
        'Color',colors(k,:), 'Marker',markers.(sev), ...
        'DisplayName',sev);
end
hold off;
xlabel('Time Lag [s]'); ylabel('Normalized ACF [rad^2]');
title('Diffractive Phase Residuals ACF');
% legend('Location','best', 'Direction', 'reverse'); 
grid on;
set(gca, 'FontSize', font_size);

%% Export the residuals ACFs plots and CSV

% Export figure as cropped, fully‐vector PDF
fig_name = 'residuals_acf_singlefreq_cpssm';  % ← change this to something meaningful
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

%% Compute the amplitude and phases periodograms and AR model PSDs

% Parameters
nfft             = 2^16;
fs               = 1/sampling_interval;
num_realizations = 300;
N                = simulation_time * fs;
win              = hamming(N);
noverlap         = 0;

% Preallocate struct for PSD results
psd_comparison = struct( ...
    'freq',     [], ...
    'amplitude', struct('periodogram',[],'ar_psd',[]), ...
    'total_phs', struct('periodogram',[],'ar_psd',[]), ...
    'refr_phs',  struct('periodogram',[],'ar_psd',[]), ...
    'diff_phs',  struct('periodogram',[],'ar_psd',[]) ...
    );

% Monte Carlo loop
seed = 2;

for i = 1:numel(severities)
    sev = severities{i};

    % initialize accumulators
    acc_per_amp   = zeros(nfft/2+1,1);
    acc_ar_amp   = zeros(nfft/2+1,1);
    acc_per_tot   = zeros(nfft/2+1,1);
    acc_ar_tot   = zeros(nfft/2+1,1);
    acc_per_refr  = zeros(nfft/2+1,1);
    acc_ar_refr  = zeros(nfft/2+1,1);
    acc_per_diff  = zeros(nfft/2+1,1);
    acc_ar_diff  = zeros(nfft/2+1,1);

    for mc = 1:num_realizations
        rng(seed + mc);

        % Generate & center signals
        [scint_ts, refr_phs_ts] = get_tppsm_multifreq_data(cpssm_params.(sev),'seed',seed+mc);
        amp_ts       = abs(scint_ts(:,1));
        total_phs_ts = get_corrected_phase(scint_ts(:,1));
        diff_phs_ts  = wrapToPi(total_phs_ts - refr_phs_ts(:,1));

        % Get periodogram using `cpsd` function
        [per_amp, f]  = cpsd(amp_ts,   amp_ts,   win, noverlap, nfft, fs);
        [per_tot, ~]  = cpsd(total_phs_ts,   total_phs_ts,   win, noverlap, nfft, fs);
        [per_refr,~]  = cpsd(refr_phs_ts(:,1),  refr_phs_ts(:,1),  win, noverlap, nfft, fs);
        [per_diff,~]  = cpsd(diff_phs_ts,  diff_phs_ts,  win, noverlap, nfft, fs);

        if isempty(psd_comparison.freq)
            psd_comparison.freq = f;
        end

        % Get AR models PSDs
        p_amp          = orders(min_sbc_idx_amp(i));
        p_tot          = orders(min_sbc_idx_total_phs(i));
        p_refr         = orders(min_sbc_idx_refr_phs(i));
        p_diff         = orders(min_sbc_idx_diff_phs(i));
        [~, A_amp, C_amp]     = arfit(amp_ts,        p_amp, p_amp);
        [~, A_tot, C_tot]     = arfit(total_phs_ts,  p_tot, p_tot);
        [~, A_refr, C_refr]   = arfit(refr_phs_ts(:,1),   p_refr,p_refr);
        [~, A_diff,C_diff]    = arfit(diff_phs_ts,   p_diff,p_diff);

        a_amp   = [1, -reshape(A_amp,1,[])];
        a_tot   = [1, -reshape(A_tot,1,[])];
        a_refr  = [1, -reshape(A_refr,1,[])];
        a_diff  = [1, -reshape(A_diff,1,[])];

        z            = exp(1j*2*pi*f*sampling_interval);
        H_amp        = 1 ./ sum(a_amp .* z.^(-(0:p_amp)), 2);
        H_tot        = 1 ./ sum(a_tot .* z.^(-(0:p_tot)), 2);
        H_refr       = 1 ./ sum(a_refr .* z.^(-(0:p_refr)), 2);
        H_diff       = 1 ./ sum(a_diff .* z.^(-(0:p_diff)), 2);

        S_ar_amp    = 2 * (C_amp/fs)  .* abs(H_amp).^2;
        S_ar_tot    = 2 * (C_tot/fs)  .* abs(H_tot).^2;
        S_ar_refr   = 2 * (C_refr/fs) .* abs(H_refr).^2;
        S_ar_diff   = 2 * (C_diff/fs) .* abs(H_diff).^2;

        % Accumulate
        acc_per_amp  = acc_per_amp  + real(per_amp);
        acc_ar_amp  = acc_ar_amp  + S_ar_amp;
        acc_per_tot  = acc_per_tot  + real(per_tot);
        acc_ar_tot  = acc_ar_tot  + S_ar_tot;
        acc_per_refr = acc_per_refr + real(per_refr);
        acc_ar_refr = acc_ar_refr + S_ar_refr;
        acc_per_diff = acc_per_diff + real(per_diff);
        acc_ar_diff = acc_ar_diff + S_ar_diff;
    end

    % Average over Monte Carlo
    psd_comparison.amplitude.(sev).periodogram  = acc_per_amp / num_realizations;
    psd_comparison.amplitude.(sev).ar_psd      = acc_ar_amp / num_realizations;
    psd_comparison.total_phs.(sev).periodogram  = acc_per_tot / num_realizations;
    psd_comparison.total_phs.(sev).ar_psd      = acc_ar_tot / num_realizations;
    psd_comparison.refr_phs.(sev).periodogram   = acc_per_refr / num_realizations;
    psd_comparison.refr_phs.(sev).ar_psd       = acc_ar_refr / num_realizations;
    psd_comparison.diff_phs.(sev).periodogram   = acc_per_diff / num_realizations;
    psd_comparison.diff_phs.(sev).ar_psd       = acc_ar_diff / num_realizations;

    seed = seed + 1;
end

%% Plot the signals periodgrams and the AR PSDs
% Plot 2×2 subplots
cmap_p = winter(numel(severities));
cmap_v = autumn(numel(severities));
metrics = {'amplitude','total_phs','refr_phs','diff_phs'};
titles  = {'Amplitude','Total Phase','Refractive Phase','Diffractive Phase'};

figure('Position',[100,100,1000,600]);
F = psd_comparison.freq;

for m = 1:4
    subplot(2,2,m); hold on;
    for k = 1:numel(severities)
        sev = severities{k};
        plot(F, 10*log10(psd_comparison.(metrics{m}).(sev).periodogram),'--', ...
            'Color',cmap_p(k,:), 'DisplayName',[sev, ' - Periodogram']);
        plot(F, 10*log10(psd_comparison.(metrics{m}).(sev).ar_psd),      '-', ...
            'Color',cmap_v(k,:), 'LineWidth',1.5, 'DisplayName', [sev, ' - AR PSD']);
    end
    hold off;
    set(gca,'XScale','log','XLim',[1e-4*fs,0.4*fs], 'FontSize', font_size);
    xlabel('Norm. freq. (× 1/T_I) [Hz]'); 
    if m > 1
        ylabel('PSD [dB(rad^2/Hz)]');
    else
        ylabel('PSD [dB/Hz]');
    end
    
    title(titles{m});
    grid on; 
    if m == 1
        legend('Location','best');
    end
end

%% Export periodograms and AR model PSDs plots and CSV
fig_name = 'periodogram_vs_ar_psd_singlefreq_cpssm';
exportgraphics(gcf, fullfile(fig_dir, [fig_name,'.pdf']), 'ContentType','vector');

T_psd = table( ...
    F, ...
    psd_comparison.amplitude.Weak.periodogram, ...
    psd_comparison.amplitude.Moderate.periodogram, ...
    psd_comparison.amplitude.Strong.periodogram, ...
    psd_comparison.amplitude.Weak.ar_psd, ...
    psd_comparison.amplitude.Moderate.ar_psd, ...
    psd_comparison.amplitude.Strong.ar_psd, ...
    psd_comparison.total_phs.Weak.periodogram, ...
    psd_comparison.total_phs.Moderate.periodogram, ...
    psd_comparison.total_phs.Strong.periodogram, ...
    psd_comparison.total_phs.Weak.ar_psd, ...
    psd_comparison.total_phs.Moderate.ar_psd, ...
    psd_comparison.total_phs.Strong.ar_psd, ...
    psd_comparison.refr_phs.Weak.periodogram, ...
    psd_comparison.refr_phs.Moderate.periodogram, ...
    psd_comparison.refr_phs.Strong.periodogram, ...
    psd_comparison.refr_phs.Weak.ar_psd, ...
    psd_comparison.refr_phs.Moderate.ar_psd, ...
    psd_comparison.refr_phs.Strong.ar_psd, ...
    psd_comparison.diff_phs.Weak.periodogram, ...
    psd_comparison.diff_phs.Moderate.periodogram, ...
    psd_comparison.diff_phs.Strong.periodogram, ...
    psd_comparison.diff_phs.Weak.ar_psd, ...
    psd_comparison.diff_phs.Moderate.ar_psd, ...
    psd_comparison.diff_phs.Strong.ar_psd, ...
    'VariableNames', { ...
    'Freq_Hz', ...
    'Amp_Weak_Per','Amp_Mod_Per','Amp_Str_Per','Amp_Weak_AR','Amp_Mod_AR','Amp_Str_AR', ...
    'Tot_Weak_Per','Tot_Mod_Per','Tot_Str_Per','Tot_Weak_AR','Tot_Mod_AR','Tot_Str_AR', ...
    'Refr_Weak_Per','Refr_Mod_Per','Refr_Str_Per','Refr_Weak_AR','Refr_Mod_AR','Refr_Str_AR', ...
    'Diff_Weak_Per','Diff_Mod_Per','Diff_Str_Per','Diff_Weak_AR','Diff_Mod_AR','Diff_Str_AR' ...
    } );

writetable(T_psd, fullfile(csv_dir, [fig_name,'.csv']));

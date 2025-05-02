% script_cpssm_var_singlefreq_fit.m
%
% Script to identify the goodness of fit of the VAR model for 
% single-frequency scintillation ampltiude and phase time series 
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

%% Monte Carlo optimal VAR model order assessment
mc_runs    = 5;
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
        [scint_ts, refr_phs_ts_training] = get_tppsm_data(cpssm_params.(severity), 'seed', seed);
        amp_ts    = abs(scint_ts);
        % NOTE: `get_corrected_phase` is a function from the
        % 'gnss_scintillation_simulator' git submodule.
        total_phs_ts_training  = get_corrected_phase(scint_ts);
        % Get the diffractive phase as the wrapped version of the
        % difference between the propagated field's phase and the
        % refractive phase time series at the ionospheric piercing point
        % (IPP).
        % NOTE: The `wrapToPi` function serves to wrap the diffractive phase
        % time series within  the [-pi, pi] bounds, in order to be
        % comparable with the CSM model in a feasible way.
        diff_phs_ts_training = wrapToPi(total_phs_ts_training - refr_phs_ts_training);

        % Fit the scintillation amplitude, total, refractive and
        % diffractive phases, respectively to a AR model, respectively.
        [~, A_amp, ~, sbc_amp]   = arfit(amp_ts,   min_order, max_order);
        [~, A_total_phs, ~, sbc_total_phs] = arfit(total_phs_ts_training, min_order, max_order);
        [~, A_refr_phs, ~, sbc_refr_phs] = arfit(refr_phs_ts_training, min_order, max_order);
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

figure('Position',[50,50,1100,400]);
colors = lines(numel(severities));  % or get(gca,'ColorOrder')

% Amplitude
subplot(2,2,1);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_amp(:,j), '-o', ...
         'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
         'DisplayName',severities{j});
end
hold off;
xlabel('VAR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal VAR Order – Amplitude');
legend('Location','best');
grid on;

% Total phase
subplot(2,2,2);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_total_phs(:,j), '-o', ...
         'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
         'DisplayName',severities{j});
end
hold off;
xlabel('VAR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal VAR Order – Total Phase');
legend('Location','best');
grid on;

% Refractive phase
subplot(2,2,3);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_refr_phs(:,j), '-o', ...
         'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
         'DisplayName',severities{j});
end
hold off;
xlabel('VAR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal VAR Order – Refractive Phase');
legend('Location','best');
grid on;

% Diffractive phase
subplot(2,2,4);
hold on;
for j = 1:numel(severities)
    plot(orders, pct_diff_phs(:,j), '-o', ...
         'LineWidth',1.5, 'MarkerSize',6, 'Color',colors(j,:), ...
         'DisplayName',severities{j});
end
hold off;
xlabel('VAR Model Order');
ylabel('Percentage of runs [%]');
title('Optimal VAR Order – Diffractive Phase');
legend('Location','best');
grid on;

% Export Optimal VAR Order Plot and CSV -----------------------------------

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

% Plot Mean SBC for Amplitude & Phase Components --------------------------
orders               = min_order:max_order;
mean_sbc_amp         = squeeze(mean(sbc_amp_array,         1)).';   % [severity × order]
mean_sbc_total_phs   = squeeze(mean(sbc_total_phs_array,   1)).';
mean_sbc_refr_phs    = squeeze(mean(sbc_refr_phs_array,    1)).';
mean_sbc_diff_phs    = squeeze(mean(sbc_diff_phs_array,    1)).';

figure('Position',[50,50,1100,400]);

% Amplitude
subplot(2,2,1);
colors = get(gca,'ColorOrder');
plot(orders, mean_sbc_amp, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_amp,2)
    [minVal, minIdx] = min(mean_sbc_amp(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('VAR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Amplitude');
legend(severities,'Location','best'); grid on; hold off;

% Total phase
subplot(2,2,2);
plot(orders, mean_sbc_total_phs, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_total_phs,2)
    [minVal, minIdx] = min(mean_sbc_total_phs(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('VAR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Total Phase');
legend(severities,'Location','best'); grid on; hold off;

% Refractive phase
subplot(2,2,3);
plot(orders, mean_sbc_refr_phs, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_refr_phs,2)
    [minVal, minIdx] = min(mean_sbc_refr_phs(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('VAR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Refractive Phase');
legend(severities,'Location','best'); grid on; hold off;

% Diffractive phase
subplot(2,2,4);
plot(orders, mean_sbc_diff_phs, 'LineWidth',1.5); hold on;
for j = 1:size(mean_sbc_diff_phs,2)
    [minVal, minIdx] = min(mean_sbc_diff_phs(:,j));
    plot(orders(minIdx), minVal, '*', 'MarkerSize',10, 'Color',colors(j,:));
end
xlabel('VAR Model Order'); ylabel('Mean SBC');
title('Mean SBC – Diffractive Phase');
legend(severities,'Location','best'); grid on; hold off;

%% Export Mean SBC Plot and CSV

% Export Mean SBC figure as a cropped, fully-vector PDF
fig_name = 'mean_sbc';
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

%% Obtain residuals for most frequent orders
[~, highest_freq_idx_amp]   = max(counts_amp,[],1);
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
    [scint_ts, refr_phs_ts_training] = get_tppsm_data(cpssm_params.(severity), 'seed', seed);
    amp_ts_training    = abs(scint_ts);
    total_phs_ts_training  = get_corrected_phase(scint_ts);
    diff_phs_ts_training = wrapToPi(total_phs_ts_training - refr_phs_ts_training);

    ord_amp = orders(highest_freq_idx_amp(i));
    ord_total_phs = orders(highest_freq_idx_total_phs(i));
    ord_refr_phs = orders(highest_freq_idx_refr_phs(i));
    ord_diff_phase = orders(highest_freq_idx_diff_phs(i));

    [w_amp, A_amp]   = arfit(amp_ts_training, ord_amp, ord_amp);
    [w_total_phs, A_total_phs] = arfit(total_phs_ts_training, ord_total_phs, ord_total_phs);
    [w_refr_phs, A_refr_phs] = arfit(refr_phs_ts_training, ord_refr_phs, ord_refr_phs);
    [w_diff_phase, A_diff_phase] = arfit(diff_phs_ts_training, ord_diff_phase, ord_diff_phase);

    % Obtaining the residuals
    % NOTE: Another generation seed is used herein to evaluate the reisuduals 
    % correlation with statistical significance.
    rng(seed+1);
    [scint_ts, refr_phs_ts_training] = get_tppsm_data(cpssm_params.(severity), 'seed', seed + 1);
    amp_ts_for_residue    = abs(scint_ts);
    total_phs_ts_for_residue  = get_corrected_phase(scint_ts);
    diff_phs_ts_for_residue = wrapToPi(total_phs_ts_for_residue - refr_phs_ts_training);

    [~, res_amp]   = arres(w_amp,   A_amp,   amp_ts_for_residue, 60);
    [~, res_total_phs]   = arres(w_total_phs,   A_total_phs,   total_phs_ts_for_residue, 60);
    [~, res_refr_phs]   = arres(w_refr_phs,   A_refr_phs,   refr_phs_ts_training, 60);
    [~, res_diff_phase]   = arres(w_diff_phase,   A_diff_phase,   diff_phs_ts_for_residue, 60);

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

% Plot the downsampled residuals of the fitted model ----------------------------------

% Define plot order and time axis
plot_order = {'Strong','Moderate','Weak'};

% Invert default colormap so
% strong→yellow (3rd), moderate→red (2nd), weak→blue (1st)
base_colors = lines(numel(plot_order));
colors      = base_colors([3,2,1],:);

figure('Position',[50,50,1100,600]);

% 1) Amplitude Residuals
subplot(2,2,1); hold on;
for k = 1:3
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.amplitude.(sev), ...
         'LineWidth',1, 'Color',colors(k,:), ...
         'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals');
title('Amplitude Residuals');
legend({'Strong','Moderate','Weak'},'Location','best');
grid on;

% 2) Total Phase Residuals
subplot(2,2,2); hold on;
for k = 1:3
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.total_phs.(sev), ...
         'LineWidth',1, 'Color',colors(k,:), ...
         'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals');
title('Total Phase Residuals');
legend({'Strong','Moderate','Weak'},'Location','best');
grid on;

% 3) Refractive Phase Residuals
subplot(2,2,3); hold on;
for k = 1:3
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.refr_phs.(sev), ...
         'LineWidth',1, 'Color',colors(k,:), ...
         'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals');
title('Refractive Phase Residuals');
legend({'Strong','Moderate','Weak'},'Location','best');
grid on;

% 4) Diffractive Phase Residuals
subplot(2,2,4); hold on;
for k = 1:3
    sev = plot_order{k};
    plot(time_downsamp, residuals_downsamp.diff_phs.(sev), ...
         'LineWidth',1, 'Color',colors(k,:), ...
         'DisplayName', sev);
end
hold off;
xlabel('Time [s]'); ylabel('Residuals');
title('Diffractive Phase Residuals');
legend({'Strong','Moderate','Weak'},'Location','best');
grid on;

% Export Residuals Plot and CSV -------------------------------------------

% Export residuals figure as a cropped, fully-vector PDF
fig_name = 'residuals';
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

% Compute one-sided ACFs of residuals -------------------------------------

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

plot_order = {'Strong','Moderate','Weak'};
time_lag   = (0:lags) * sampling_interval;
colors     = lines(numel(plot_order));
colors     = colors([3,2,1],:);    % strong→yellow, moderate→red, weak→blue
stem_width = 1.5;
markers    = struct('Weak','o','Moderate','s','Strong','^');

% Plot ACFs in 2×2 grid: Amplitude, Total, Refractive, Diffractive
figure('Position',[50,50,1100,600]);

% 1) Amplitude ACF
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

% 2) Total Phase ACF
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

% 3) Refractive Phase ACF
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

% 4) Diffractive Phase ACF
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

%% Monte Carlo CPSD-based PSD vs. VAR PSD Comparison (with Refractive Phase)

% ─── Parameters ──────────────────────────────────────────────────────
nfft             = 2^16;
fs               = 1/sampling_interval;
num_realizations = 1;
N                = simulation_time * fs;
win              = hamming(N);
noverlap         = 0;

% ─── Preallocate struct for PSD results ──────────────────────────────
psd_comparison = struct( ...
  'freq',     [], ...
  'amplitude', struct('periodogram',[],'var_psd',[]), ...
  'total_phs', struct('periodogram',[],'var_psd',[]), ...
  'refr_phs',  struct('periodogram',[],'var_psd',[]), ...
  'diff_phs',  struct('periodogram',[],'var_psd',[]) ...
);

% ─── Monte Carlo loop ─────────────────────────────────────────────────
seed = 2;
for i = 1:numel(severities)
    sev = severities{i};

    % initialize accumulators
    acc_per_amp   = zeros(nfft/2+1,1);
    acc_var_amp   = zeros(nfft/2+1,1);
    acc_per_tot   = zeros(nfft/2+1,1);
    acc_var_tot   = zeros(nfft/2+1,1);
    acc_per_refr  = zeros(nfft/2+1,1);
    acc_var_refr  = zeros(nfft/2+1,1);
    acc_per_diff  = zeros(nfft/2+1,1);
    acc_var_diff  = zeros(nfft/2+1,1);

    for mc = 1:num_realizations
        rng(seed + mc);

        % 1) Generate & center signals
        [scint_ts, refr_phs_ts] = get_tppsm_data(cpssm_params.(sev),'seed',seed+mc);
        amp_ts       = abs(scint_ts);
        total_phs_ts = get_corrected_phase(scint_ts);
        diff_phs_ts  = wrapToPi(total_phs_ts - refr_phs_ts);

        amp_ctr   = amp_ts   - mean(amp_ts);
        tot_ctr   = total_phs_ts - mean(total_phs_ts);
        refr_ctr  = refr_phs_ts  - mean(refr_phs_ts);
        diff_ctr  = diff_phs_ts  - mean(diff_phs_ts);

        % 2) Get periodogram using `cpsd` function
        [per_amp, f]  = cpsd(amp_ctr,   amp_ctr,   win, noverlap, nfft, fs);
        [per_tot, ~]  = cpsd(tot_ctr,   tot_ctr,   win, noverlap, nfft, fs);
        [per_refr,~]  = cpsd(refr_ctr,  refr_ctr,  win, noverlap, nfft, fs);
        [per_diff,~]  = cpsd(diff_ctr,  diff_ctr,  win, noverlap, nfft, fs);

        if isempty(psd_comparison.freq)
            psd_comparison.freq = f;
        end

        % 3) Get AR models PSDs
        p_amp          = orders(highest_freq_idx_amp(i));
        p_tot          = orders(highest_freq_idx_total_phs(i));
        p_refr         = orders(highest_freq_idx_refr_phs(i));
        p_diff         = orders(highest_freq_idx_diff_phs(i));
        [~, A_amp, C_amp]     = arfit(amp_ts,        p_amp, p_amp);
        [~, A_tot, C_tot]     = arfit(total_phs_ts,  p_tot, p_tot);
        [~, A_refr, C_refr]   = arfit(refr_phs_ts,   p_refr,p_refr);
        [~, A_diff,C_diff]    = arfit(diff_phs_ts,   p_diff,p_diff);

        a_amp   = [1, -reshape(A_amp,1,[])];
        a_tot   = [1, -reshape(A_tot,1,[])];
        a_refr  = [1, -reshape(A_refr,1,[])];
        a_diff  = [1, -reshape(A_diff,1,[])];

        z            = exp(-1j*2*pi*f*sampling_interval);
        H_amp        = 1 ./ sum(a_amp .* z.^(-(0:p_amp)), 2);
        H_tot        = 1 ./ sum(a_tot .* z.^(-(0:p_tot)), 2);
        H_refr       = 1 ./ sum(a_refr .* z.^(-(0:p_refr)), 2);
        H_diff       = 1 ./ sum(a_diff .* z.^(-(0:p_diff)), 2);

        S_var_amp    = 2 * (C_amp/fs)  .* abs(H_amp).^2;
        S_var_tot    = 2 * (C_tot/fs)  .* abs(H_tot).^2;
        S_var_refr   = 2 * (C_refr/fs) .* abs(H_refr).^2;
        S_var_diff   = 2 * (C_diff/fs) .* abs(H_diff).^2;

        % accumulate
        acc_per_amp  = acc_per_amp  + real(per_amp);
        acc_var_amp  = acc_var_amp  + S_var_amp;
        acc_per_tot  = acc_per_tot  + real(per_tot);
        acc_var_tot  = acc_var_tot  + S_var_tot;
        acc_per_refr = acc_per_refr + real(per_refr);
        acc_var_refr = acc_var_refr + S_var_refr;
        acc_per_diff = acc_per_diff + real(per_diff);
        acc_var_diff = acc_var_diff + S_var_diff;
    end

    % average over Monte Carlo
    psd_comparison.amplitude.(sev).periodogram  = acc_per_amp / num_realizations;
    psd_comparison.amplitude.(sev).var_psd      = acc_var_amp / num_realizations;
    psd_comparison.total_phs.(sev).periodogram  = acc_per_tot / num_realizations;
    psd_comparison.total_phs.(sev).var_psd      = acc_var_tot / num_realizations;
    psd_comparison.refr_phs.(sev).periodogram   = acc_per_refr / num_realizations;
    psd_comparison.refr_phs.(sev).var_psd       = acc_var_refr / num_realizations;
    psd_comparison.diff_phs.(sev).periodogram   = acc_per_diff / num_realizations;
    psd_comparison.diff_phs.(sev).var_psd       = acc_var_diff / num_realizations;

    seed = seed + 1;
end

% ─── Plot 2×2 Subplots ──────────────────────────────────────────────────
cmap_p = winter(numel(severities));
cmap_v = autumn(numel(severities));
metrics = {'amplitude','total_phs','refr_phs','diff_phs'};
titles  = {'Amplitude','Total Phase','Refractive Phase','Diffractive Phase'};

figure('Position',[50,50,1200,800]);
F = psd_comparison.freq;

for m = 1:4
    subplot(2,2,m); hold on;
    for k = 1:numel(severities)
        sev = severities{k};
        plot(F, 10*log10(psd_comparison.(metrics{m}).(sev).periodogram),'--', ...
             'Color',cmap_p(k,:), 'LineWidth',1);
        plot(F, 10*log10(psd_comparison.(metrics{m}).(sev).var_psd),      '-', ...
             'Color',cmap_v(k,:), 'LineWidth',2, 'DisplayName',sev);
    end
    hold off;
    set(gca,'XScale','log','XLim',[1e-4*fs,0.4*fs]);
    xlabel('Frequency [Hz]'); ylabel('Power [dB]');
    title([titles{m} ' PSD vs. VAR PSD']);
    grid on; legend('Location','best');
end

% ─── Export Figure & CSV ───────────────────────────────────────────────
fig_name = 'psd_var_comparison';
exportgraphics(gcf, fullfile(fig_dir, [fig_name,'.pdf']), 'ContentType','vector');

T_psd = table( ...
    F, ...
    psd_comparison.amplitude.Weak.periodogram, ...
    psd_comparison.amplitude.Moderate.periodogram, ...
    psd_comparison.amplitude.Strong.periodogram, ...
    psd_comparison.amplitude.Weak.var_psd, ...
    psd_comparison.amplitude.Moderate.var_psd, ...
    psd_comparison.amplitude.Strong.var_psd, ...
    psd_comparison.total_phs.Weak.periodogram, ...
    psd_comparison.total_phs.Moderate.periodogram, ...
    psd_comparison.total_phs.Strong.periodogram, ...
    psd_comparison.total_phs.Weak.var_psd, ...
    psd_comparison.total_phs.Moderate.var_psd, ...
    psd_comparison.total_phs.Strong.var_psd, ...
    psd_comparison.refr_phs.Weak.periodogram, ...
    psd_comparison.refr_phs.Moderate.periodogram, ...
    psd_comparison.refr_phs.Strong.periodogram, ...
    psd_comparison.refr_phs.Weak.var_psd, ...
    psd_comparison.refr_phs.Moderate.var_psd, ...
    psd_comparison.refr_phs.Strong.var_psd, ...
    psd_comparison.diff_phs.Weak.periodogram, ...
    psd_comparison.diff_phs.Moderate.periodogram, ...
    psd_comparison.diff_phs.Strong.periodogram, ...
    psd_comparison.diff_phs.Weak.var_psd, ...
    psd_comparison.diff_phs.Moderate.var_psd, ...
    psd_comparison.diff_phs.Strong.var_psd, ...
    'VariableNames', { ...
      'Freq_Hz', ...
      'Amp_Weak_Per','Amp_Mod_Per','Amp_Str_Per','Amp_Weak_VAR','Amp_Mod_VAR','Amp_Str_VAR', ...
      'Tot_Weak_Per','Tot_Mod_Per','Tot_Str_Per','Tot_Weak_VAR','Tot_Mod_VAR','Tot_Str_VAR', ...
      'Refr_Weak_Per','Refr_Mod_Per','Refr_Str_Per','Refr_Weak_VAR','Refr_Mod_VAR','Refr_Str_VAR', ...
      'Diff_Weak_Per','Diff_Mod_Per','Diff_Str_Per','Diff_Weak_VAR','Diff_Mod_VAR','Diff_Str_VAR' ...
    } );

writetable(T_psd, fullfile(csv_dir, [fig_name,'.csv']));

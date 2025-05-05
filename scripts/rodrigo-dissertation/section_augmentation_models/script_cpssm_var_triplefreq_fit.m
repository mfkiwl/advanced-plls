%% script_cpssm_var_triplefreq_fit_separate.m
% Separate VAR modeling for Triple-Frequency Amplitudes and Total Phases
% using the CPSSM and ARFIT algorithm [1].
%
% References:
% [1] Schneider, Tapio, and Arnold Neumaier. “Algorithm 808: ARfit—a Matlab
% Package for the Estimation of Parameters and Eigenmodes of Multivariate
% Autoregressive Models.” ACM Trans. Math. Softw. 27, no. 1 (2001): 58–65.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;
addpath(genpath(fullfile(pwd,'..','..','..','libs')));

% Output folders
fig_dir = 'pdf_figures_cpssm_triplefreq';
csv_dir = 'csv_data_cpssm_triplefreq';
if ~exist(fig_dir,'dir'),  mkdir(fig_dir);  end
if ~exist(csv_dir,'dir'),  mkdir(csv_dir);  end

%% Simulation parameters
simulation_time    = 300;
sampling_interval  = 0.01;
severities         = {'Weak','Moderate','Strong'};
frequency_bands    = {'L1','L2','L5'};
freq_amount        = numel(frequency_bands);
cpssm_params = struct( ...
    'Weak',    {'weak',     'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 1.5},...
    'Moderate',{'moderate', 'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 0.8},...
    'Strong',  {'strong',   'is_enable_cmd_print', false, 'simulation_time', simulation_time, 'sampling_interval', sampling_interval, 'rhof_veff_ratio', 0.27}...
    );
%% Monte Carlo optimal order assessment
mc_runs    = 10;
min_order  = 1;
max_order  = 30;
orders_vec = min_order:max_order;
n_orders   = numel(orders_vec);

% Preallocate
optimal_orders_amp = zeros(mc_runs, numel(severities));
sbc_amp_array      = zeros(mc_runs, numel(severities), n_orders);
optimal_orders_phs = zeros(mc_runs, numel(severities));
sbc_phs_array      = zeros(mc_runs, numel(severities), n_orders);

seed = 1;
for mc_idx = 1:mc_runs
    for i = 1:numel(severities)
        sev = severities{i};
        rng(seed);
        [scint_ts, ~] = get_tppsm_data(cpssm_params.(sev), 'seed', seed);
        amp_ts = abs(scint_ts);
        phs_ts = zeros(size(amp_ts));
        for f = 1:freq_amount
            phs_ts(:,f) = get_corrected_phase(scint_ts(:,f));
        end
        % Fit VAR to amplitudes
        [~, A_amp, ~, sbc_amp] = arfit(amp_ts, min_order, max_order);
        optimal_orders_amp(mc_idx,i) = size(A_amp,2)/size(A_amp,1);
        sbc_amp_array(mc_idx,i,:)      = sbc_amp;
        % Fit VAR to total phases
        [~, A_phs, ~, sbc_phs] = arfit(phs_ts, min_order, max_order);
        optimal_orders_phs(mc_idx,i) = size(A_phs,2)/size(A_phs,1);
        sbc_phs_array(mc_idx,i,:)      = sbc_phs;
        seed = seed + 1;
    end
end

% Histogram counts
counts_amp = zeros(n_orders, numel(severities));
counts_phs = zeros(n_orders, numel(severities));
for i = 1:numel(severities)
    counts_amp(:,i) = histcounts(optimal_orders_amp(:,i), [orders_vec, orders_vec(end)+1]);
    counts_phs(:,i) = histcounts(optimal_orders_phs(:,i), [orders_vec, orders_vec(end)+1]);
end
% Convert to percent
pct_amp = counts_amp / mc_runs * 100;
pct_phs = counts_phs / mc_runs * 100;

%% Plot optimal order frequencies
colors = lines(numel(severities));

% Amplitude orders
figure('Position',[50,50,800,600]);
subplot(2,1,1);
hold on;
for j = 1:numel(severities)
    plot(orders_vec, pct_amp(:,j), '-o', 'LineWidth',1.5, 'MarkerSize',6, ...
         'Color',colors(j,:), 'DisplayName',severities{j});
end
hold off;
xlabel('VAR Order'); ylabel('Percentage [%]');
title('Optimal VAR Order – Amplitudes');
legend('Location','best'); grid on;
exportgraphics(gcf, fullfile(fig_dir,'opt_order_amp.pdf'),'ContentType','vector');
% CSV
T_amp = table(orders_vec.', pct_amp(:,1), pct_amp(:,2), pct_amp(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_amp, fullfile(csv_dir,'opt_order_amp.csv'));

% Phase orders
subplot(2,1,2);
hold on;
for j = 1:numel(severities)
    plot(orders_vec, pct_phs(:,j), '-s', 'LineWidth',1.5, 'MarkerSize',6, ...
         'Color',colors(j,:), 'DisplayName',severities{j});
end
hold off;
xlabel('VAR Order'); ylabel('Percentage [%]');
title('Optimal VAR Order – Total Phases');
legend('Location','best'); grid on;
exportgraphics(gcf, fullfile(fig_dir,'opt_order_phs.pdf'),'ContentType','vector');
% CSV
T_phs = table(orders_vec.', pct_phs(:,1), pct_phs(:,2), pct_phs(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_phs, fullfile(csv_dir,'opt_order_phs.csv'));

%% Compute mean SBC and plot
mean_sbc_amp = squeeze(mean(sbc_amp_array,1)).';
mean_sbc_phs = squeeze(mean(sbc_phs_array,1)).';

% Amplitude SBC
figure('Position',[50,50,800,600]);
subplot(2,1,1);
hold on;
plot(orders_vec, mean_sbc_amp, 'LineWidth',1.5);
for j=1:size(mean_sbc_amp,2)
    [minv, idx] = min(mean_sbc_amp(:,j));
    plot(orders_vec(idx),minv,'*','MarkerSize',10,'Color',colors(j,:));
end
hold off;
xlabel('VAR Order'); ylabel('Mean SBC');
title('Mean SBC – Amplitudes');
legend(severities,'Location','best'); grid on;
exportgraphics(gcf, fullfile(fig_dir,'mean_sbc_amp.pdf'),'ContentType','vector');
% CSV
T_sbc_amp = table(orders_vec.', mean_sbc_amp(:,1), mean_sbc_amp(:,2), mean_sbc_amp(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_sbc_amp, fullfile(csv_dir,'mean_sbc_amp.csv'));

% Phase SBC
subplot(2,1,2);
hold on;
plot(orders_vec, mean_sbc_phs, 'LineWidth',1.5);
for j=1:size(mean_sbc_phs,2)
    [minv, idx] = min(mean_sbc_phs(:,j));
    plot(orders_vec(idx),minv,'*','MarkerSize',10,'Color',colors(j,:));
end
hold off;
xlabel('VAR Order'); ylabel('Mean SBC');
title('Mean SBC – Total Phases');
legend(severities,'Location','best'); grid on;
exportgraphics(gcf, fullfile(fig_dir,'mean_sbc_phs.pdf'),'ContentType','vector');
% CSV
T_sbc_phs = table(orders_vec.', mean_sbc_phs(:,1), mean_sbc_phs(:,2), mean_sbc_phs(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_sbc_phs, fullfile(csv_dir,'mean_sbc_phs.csv'));

%% Residual analysis for separate VAR fits
[~, idx_amp] = max(counts_amp,[],1);
[~, idx_phs] = max(counts_phs,[],1);

residuals = struct('amplitude',[],'total_phase',[]);
residuals_down = struct('amplitude',[],'total_phase',[]);

seed = 100;
for i = 1:numel(severities)
    sev = severities{i};
    rng(seed);
    [scint_train,~] = get_tppsm_data(cpssm_params.(sev),'seed',seed);
    amp_train = abs(scint_train);
    phs_train = zeros(size(amp_train));
    for f=1:freq_amount, phs_train(:,f)=get_corrected_phase(scint_train(:,f)); end
    ord_a = orders_vec(idx_amp(i));
    ord_p = orders_vec(idx_phs(i));
    [w_a,A_a] = arfit(amp_train, ord_a, ord_a);
    [w_p,A_p] = arfit(phs_train, ord_p, ord_p);
    % generate fresh data for residuals
    rng(seed+1);
    [scint_res,~] = get_tppsm_data(cpssm_params.(sev),'seed',seed+1);
    amp_res_ts = abs(scint_res);
    phs_res_ts = zeros(size(amp_res_ts)); for f=1:freq_amount, phs_res_ts(:,f)=get_corrected_phase(scint_res(:,f)); end
    [~, res_a] = arres(w_a,A_a, amp_res_ts, 60);
    [~, res_p] = arres(w_p,A_p, phs_res_ts, 60);
    % pad NaNs
    residuals.amplitude.(sev)    = [NaN(ord_a,size(res_a,2)); res_a];
    residuals.total_phase.(sev) = [NaN(ord_p,size(res_p,2)); res_p];
    seed = seed + 1;
end

%% Downsample, plot residuals (3×2 layout), and export CSVs
time_vec = sampling_interval : sampling_interval : simulation_time;
down     = 1;
time_ds  = downsample(time_vec, down);

for i = 1:numel(severities)
    sev = severities{i};
    residuals_down.amplitude.(sev)    = downsample(residuals.amplitude.(sev), down);
    residuals_down.total_phase.(sev) = downsample(residuals.total_phase.(sev), down);
end

plot_order = {'Strong','Moderate','Weak'};
bands      = {'L1','L2','L5'};
colors     = lines(numel(plot_order));

figure('Position',[100,100,1200,900]);
for r = 1:3
    % Amplitude residuals
    subplot(3,2,(r-1)*2+1); hold on;
    for k = 1:3
        sev = plot_order{k};
        plot(time_ds, residuals_down.amplitude.(sev)(:,r), ...
             'LineWidth',1, 'Color',colors(k,:), 'DisplayName',sev);
    end
    hold off;
    xlabel('Time [s]'); ylabel('Residual');
    title(sprintf('Amplitude Residuals (%s)', bands{r}));
    legend('Location','best'); grid on;

    % Total phase residuals
    subplot(3,2,(r-1)*2+2); hold on;
    for k = 1:3
        sev = plot_order{k};
        plot(time_ds, residuals_down.total_phase.(sev)(:,r), ...
             'LineWidth',1, 'Color',colors(k,:), 'DisplayName',sev);
    end
    hold off;
    xlabel('Time [s]'); ylabel('Residual');
    title(sprintf('Total Phase Residuals (%s)', bands{r}));
    legend('Location','best'); grid on;
end

% Export figure
exportgraphics(gcf, fullfile(fig_dir,'residuals_separate.pdf'),'ContentType','vector');

% Build and export CSV of residuals
% Columns: Time_s, Amp_L1_Weak, Amp_L1_Moderate, Amp_L1_Strong, ..., Tot_L5_Strong
T_res = table( ...
    time_ds.', ...
    residuals_down.amplitude.Weak(:,1),    residuals_down.amplitude.Moderate(:,1),    residuals_down.amplitude.Strong(:,1),  ...
    residuals_down.total_phase.Weak(:,1),  residuals_down.total_phase.Moderate(:,1),  residuals_down.total_phase.Strong(:,1),...
    residuals_down.amplitude.Weak(:,2),    residuals_down.amplitude.Moderate(:,2),    residuals_down.amplitude.Strong(:,2),  ...
    residuals_down.total_phase.Weak(:,2),  residuals_down.total_phase.Moderate(:,2),  residuals_down.total_phase.Strong(:,2),...
    residuals_down.amplitude.Weak(:,3),    residuals_down.amplitude.Moderate(:,3),    residuals_down.amplitude.Strong(:,3),  ...
    residuals_down.total_phase.Weak(:,3),  residuals_down.total_phase.Moderate(:,3),  residuals_down.total_phase.Strong(:,3),...
    'VariableNames',{ ...
      'Time_s', ...
      'Amp_L1_Weak','Amp_L1_Moderate','Amp_L1_Strong', 'Tot_L1_Weak','Tot_L1_Moderate','Tot_L1_Strong', ...
      'Amp_L2_Weak','Amp_L2_Moderate','Amp_L2_Strong', 'Tot_L2_Weak','Tot_L2_Moderate','Tot_L2_Strong', ...
      'Amp_L5_Weak','Amp_L5_Moderate','Amp_L5_Strong', 'Tot_L5_Weak','Tot_L5_Moderate','Tot_L5_Strong' } );
writetable(T_res, fullfile(csv_dir,'residuals_separate.csv'));



%% Compute, plot one‐sided ACFs (3×2 layout), and export CSVs
lags      = 20;
acf_struct = struct('amplitude',[],'total_phase',[]);

% Compute ACFs per column
for i = 1:numel(severities)
    sev = severities{i};
    % Amplitude ACF
    mat_a = residuals.amplitude.(sev);
    acf_a = zeros(lags+1, freq_amount);
    for c = 1:freq_amount
        d = mat_a(~isnan(mat_a(:,c)), c);
        tmp = xcorr(d, d, lags, 'normalized');
        acf_a(:, c) = tmp(lags+1:end);
    end
    acf_struct.amplitude.(sev) = acf_a;
    % Total phase ACF
    mat_p = residuals.total_phase.(sev);
    acf_p = zeros(lags+1, freq_amount);
    for c = 1:freq_amount
        d = mat_p(~isnan(mat_p(:,c)), c);
        tmp = xcorr(d, d, lags, 'normalized');
        acf_p(:, c) = tmp(lags+1:end);
    end
    acf_struct.total_phase.(sev) = acf_p;
end

% Plot ACFs
stem_w  = 1.5;
markers = struct('Weak','o','Moderate','s','Strong','^');
time_lag = (0:lags) * sampling_interval;

figure('Position',[100,100,1200,900]);
for r = 1:3
    % Amplitude ACF
    subplot(3,2,(r-1)*2+1); hold on;
    for k = 1:3
        sev = plot_order{k};
        stem(time_lag, acf_struct.amplitude.(sev)(:,r), ...
             'LineWidth',stem_w, 'Marker',markers.(sev), 'DisplayName',sev);
    end
    hold off;
    xlabel('Lag [s]'); ylabel('ACF');
    title(sprintf('Amplitude ACF (%s)', bands{r}));
    legend('Location','best'); grid on;

    % Total phase ACF
    subplot(3,2,(r-1)*2+2); hold on;
    for k = 1:3
        sev = plot_order{k};
        stem(time_lag, acf_struct.total_phase.(sev)(:,r), ...
             'LineWidth',stem_w, 'Marker',markers.(sev), 'DisplayName',sev);
    end
    hold off;
    xlabel('Lag [s]'); ylabel('ACF');
    title(sprintf('Total Phase ACF (%s)', bands{r}));
    legend('Location','best'); grid on;
end

% Export figure
exportgraphics(gcf, fullfile(fig_dir,'residuals_acf_separate.pdf'),'ContentType','vector');

% Build and export CSV of ACFs
% Columns: Lag_s, Amp_L1_Strong, Amp_L1_Moderate, Amp_L1_Weak, Tot_L1_Strong, …
T_acf = table( ...
    time_lag.', ...
    acf_struct.amplitude.Strong(:,1),   acf_struct.amplitude.Moderate(:,1),   acf_struct.amplitude.Weak(:,1),   ...
    acf_struct.total_phase.Strong(:,1), acf_struct.total_phase.Moderate(:,1), acf_struct.total_phase.Weak(:,1), ...
    acf_struct.amplitude.Strong(:,2),   acf_struct.amplitude.Moderate(:,2),   acf_struct.amplitude.Weak(:,2),   ...
    acf_struct.total_phase.Strong(:,2), acf_struct.total_phase.Moderate(:,2), acf_struct.total_phase.Weak(:,2), ...
    acf_struct.amplitude.Strong(:,3),   acf_struct.amplitude.Moderate(:,3),   acf_struct.amplitude.Weak(:,3),   ...
    acf_struct.total_phase.Strong(:,3), acf_struct.total_phase.Moderate(:,3), acf_struct.total_phase.Weak(:,3), ...
    'VariableNames',{ ...
      'Lag_s', ...
      'Amp_L1_Strong','Amp_L1_Moderate','Amp_L1_Weak', 'Tot_L1_Strong','Tot_L1_Moderate','Tot_L1_Weak', ...
      'Amp_L2_Strong','Amp_L2_Moderate','Amp_L2_Weak', 'Tot_L2_Strong','Tot_L2_Moderate','Tot_L2_Weak', ...
      'Amp_L5_Strong','Amp_L5_Moderate','Amp_L5_Weak', 'Tot_L5_Strong','Tot_L5_Moderate','Tot_L5_Weak' } );
writetable(T_acf, fullfile(csv_dir,'residuals_acf_separate.csv'));

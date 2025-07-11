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

font_size = 11;

%% Monte Carlo optimal order assessment
mc_runs    = 300;
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
        [scint_ts, ~] = get_tppsm_multifreq_data(cpssm_params.(sev), 'seed', seed);
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
figure('Position',[100,100,1000,250]);
subplot(1,2,1);
hold on;
for j = 1:numel(severities)
    plot(orders_vec, pct_amp(:,j), '-o', 'LineWidth',1.5, 'MarkerSize',6, ...
         'Color',colors(j,:), 'DisplayName',severities{j});
end
hold off;
xlabel('VAR Order'); ylabel('Percentage [%]');
title('Optimal VAR Order – Amplitudes');
legend('Location','best'); grid on;
set(gca, 'FontSize', font_size);
% CSV
T_amp = table(orders_vec.', pct_amp(:,1), pct_amp(:,2), pct_amp(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_amp, fullfile(csv_dir,'optimal_var_order_frequency_amp_triplefreq_cpssm.csv'));

% Phase orders
subplot(1,2,2);
hold on;
for j = 1:numel(severities)
    plot(orders_vec, pct_phs(:,j), '-s', 'LineWidth',1.5, 'MarkerSize',6, ...
         'Color',colors(j,:), 'DisplayName',severities{j});
end
hold off;
xlabel('VAR Order'); ylabel('Percentage [%]');
title('Optimal VAR Order – Total Phases');
% legend('Location','best');
grid on;
set(gca, 'FontSize', font_size);

exportgraphics(gcf, fullfile(fig_dir,'optimal_var_order_frequency_triplefreq_cpssm.pdf'),'ContentType','vector');
% CSV
T_phs = table(orders_vec.', pct_phs(:,1), pct_phs(:,2), pct_phs(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_phs, fullfile(csv_dir,'opt_order.csv'));

%% Compute mean SBC and plot
mean_sbc_amp = squeeze(mean(sbc_amp_array,1)).';
mean_sbc_phs = squeeze(mean(sbc_phs_array,1)).';

% Amplitude SBC
figure('Position',[100,100,1000,250]);
subplot(1,2,1);
hold on;
plot(orders_vec, mean_sbc_amp, 'LineWidth',1.5);
for j=1:size(mean_sbc_amp,2)
    [minv, idx] = min(mean_sbc_amp(:,j));
    plot(orders_vec(idx),minv,'*','MarkerSize',10,'Color',colors(j,:));
end
hold off;
xlabel('VAR Order'); ylabel('Mean SBC');
title('Mean SBC – Amplitudes');
legend(severities,'Location','northeast'); grid on;
set(gca, 'FontSize', font_size);

%exportgraphics(gcf, fullfile(fig_dir,'mean_sbc_amp.pdf'),'ContentType','vector');
% CSV
T_sbc_amp = table(orders_vec.', mean_sbc_amp(:,1), mean_sbc_amp(:,2), mean_sbc_amp(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_sbc_amp, fullfile(csv_dir,'mean_sbc_amp_triplefreq_cpssm.csv'));

% Phase SBC
subplot(1,2,2);
hold on;
plot(orders_vec, mean_sbc_phs, 'LineWidth',1.5);
for j=1:size(mean_sbc_phs,2)
    [minv, idx] = min(mean_sbc_phs(:,j));
    plot(orders_vec(idx),minv,'*','MarkerSize',10,'Color',colors(j,:));
end
hold off;
xlabel('VAR Order'); ylabel('Mean SBC');
title('Mean SBC – Total Phases');
% legend(severities,'Location','best');
grid on;
set(gca, 'FontSize', font_size);

exportgraphics(gcf, fullfile(fig_dir,'mean_sbc_triplefreq_cpssm.pdf'),'ContentType','vector');
% CSV
T_sbc_phs = table(orders_vec.', mean_sbc_phs(:,1), mean_sbc_phs(:,2), mean_sbc_phs(:,3), ...
    'VariableNames',{'Order','Weak','Moderate','Strong'});
writetable(T_sbc_phs, fullfile(csv_dir,'mean_sbc_phs.csv'));

%% Residual analysis for separate VAR fits
[~, min_sbc_idx_amp] = min(mean_sbc_amp,[],1);
[~, min_sbc_idx_phs] = min(mean_sbc_phs,[],1);

residuals = struct('amplitude',[],'total_phase',[]);
residuals_down = struct('amplitude',[],'total_phase',[]);

seed = 100;
for i = 1:numel(severities)
    sev = severities{i};
    rng(seed);
    [scint_train,~] = get_tppsm_multifreq_data(cpssm_params.(sev),'seed',seed);
    amp_train = abs(scint_train);
    phs_train = zeros(size(amp_train));
    for f=1:freq_amount, phs_train(:,f)=get_corrected_phase(scint_train(:,f)); end
    ord_a = orders_vec(min_sbc_idx_amp(i));
    ord_p = orders_vec(min_sbc_idx_phs(i));
    [w_a,A_a] = arfit(amp_train, ord_a, ord_a);
    [w_p,A_p] = arfit(phs_train, ord_p, ord_p);
    % generate fresh data for residuals
    rng(seed+1);
    [scint_res,~] = get_tppsm_multifreq_data(cpssm_params.(sev),'seed',seed+1);
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
down_factor = 1;
time_ds  = downsample(time_vec, down_factor);

for i = 1:numel(severities)
    sev = severities{i};
    residuals_down.amplitude.(sev)    = downsample(residuals.amplitude.(sev), down_factor);
    residuals_down.total_phase.(sev) = downsample(residuals.total_phase.(sev), down_factor);
end

plot_order = {'Strong','Moderate','Weak'};
bands      = {'L1','L2','L5'};
colors     = flip(lines(numel(plot_order)));

figure('Position',[100,100,1000,750]);
for r = 1:3
    % Amplitude residuals
    subplot(3,2,(r-1)*2+1); hold on;
    for k = 1:3
        sev = plot_order{k};
        plot(time_ds, residuals_down.amplitude.(sev)(:,r), ...
             'LineWidth',1, 'Color',colors(k,:), 'DisplayName',sev);
    end
    hold off;
    xlabel('Time [s]');
    ylabel('Residuals');
    ylim([-0.2,0.4]);
    title(sprintf('Amplitude Residuals (%s)', bands{r}));
    if r == 2
        legend('Location','best', 'Direction','reverse');
    end
    grid on;
    set(gca, 'FontSize', font_size);
    % Total phase residuals
    subplot(3,2,(r-1)*2+2); hold on;
    for k = 1:3
        sev = plot_order{k};
        plot(time_ds, residuals_down.total_phase.(sev)(:,r), ...
             'LineWidth',1, 'Color',colors(k,:), 'DisplayName',sev);
    end
    hold off;
    xlabel('Time [s]');
    ylabel('Residuals [rad]');
    ylim([-pi,pi]);
    title(sprintf('Total Phase Residuals (%s)', bands{r}));
    % legend('Location','best', 'Direction','reverse'); grid on;
    set(gca, 'FontSize', font_size);
end

% Export figure
exportgraphics(gcf, fullfile(fig_dir,'residuals_triplefreq_cpssm.pdf'),'ContentType','vector');

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

figure('Position',[100,100,1000,750]);
for r = 1:3
    % Amplitude ACF
    subplot(3,2,(r-1)*2+1); hold on;
    for k = 1:3
        sev = plot_order{k};
        stem(time_lag, acf_struct.amplitude.(sev)(:,r), ...
             'LineWidth',stem_w, 'Marker',markers.(sev), 'DisplayName',sev, 'Color',colors(k,:));
    end
    hold off;
    xlabel('Time Lag [s]'); ylabel('Normalized ACF');
    title(sprintf('Amplitude ACF (%s)', bands{r}));
    if r == 2
        legend('Location','best', 'Direction','reverse');
    end
    grid on;
    set(gca, 'FontSize', font_size);

    % Total phase ACF
    subplot(3,2,(r-1)*2+2); hold on;
    for k = 1:3
        sev = plot_order{k};
        stem(time_lag, acf_struct.total_phase.(sev)(:,r), ...
             'LineWidth',stem_w, 'Marker',markers.(sev), 'DisplayName',sev, 'Color',colors(k,:));
    end
    hold off;
    xlabel('Time Lag [s]'); 
    ylabel('Normalized ACF [rad^2]');
    title(sprintf('Total Phase ACF (%s)', bands{r}));
    % legend('Location','best', 'Direction','reverse');
    grid on;
    set(gca, 'FontSize', font_size);
end

% Export figure
exportgraphics(gcf, fullfile(fig_dir,'residuals_acf_triplefreq_cpssm.pdf'),'ContentType','vector');

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

%% Compute auto- & cross-spectra and VAR PSDs for triple-frequency
% Parameters
nfft             = 2^16;
fs               = 1/sampling_interval;
num_realizations = 300;
N                = simulation_time * fs;
win              = hamming(N);
noverlap         = 0;

% Initialize accumulators per severity
for i = 1:numel(severities)
    sev = severities{i};
    acc_per_amp.(sev) = zeros(nfft/2+1, freq_amount, freq_amount);
    acc_var_amp.(sev) = zeros(nfft/2+1, freq_amount, freq_amount);
    acc_per_phs.(sev) = zeros(nfft/2+1, freq_amount, freq_amount);
    acc_var_phs.(sev) = zeros(nfft/2+1, freq_amount, freq_amount);
end

seed = 2;
for mc = 1:num_realizations
    for i = 1:numel(severities)
        sev = severities{i};
        rng(seed + mc);

        % Generate signals
        [scint_ts, ~] = get_tppsm_multifreq_data(cpssm_params.(sev), 'seed', seed+mc);
        amp_ts = abs(scint_ts);
        phs_ts = zeros(size(amp_ts));
        for f = 1:freq_amount
            phs_ts(:,f) = get_corrected_phase(scint_ts(:,f));
        end
        
        amp_centered_ts = amp_ts - mean(amp_ts);
        phs_centered_ts = phs_ts - mean(phs_ts);

        % auto- & cross-periodograms
        for p = 1:freq_amount
            for q = 1:freq_amount
                [P_amp, f] = cpsd(amp_ts(:,p), amp_centered_ts(:,q), win, noverlap, nfft, fs);
                acc_per_amp.(sev)(:,p,q) = acc_per_amp.(sev)(:,p,q) + abs(P_amp);
                [P_phs, ~] = cpsd(phs_ts(:,p), phs_centered_ts(:,q), win, noverlap, nfft, fs);
                acc_per_phs.(sev)(:,p,q) = acc_per_phs.(sev)(:,p,q) + abs(P_phs);
            end
        end

        % Fit VAR for amplitude & phase to get C matrices
        p_amp = orders_vec(min_sbc_idx_amp(i));
        p_phs = orders_vec(min_sbc_idx_phs(i));
        [~, A_amp, C_amp] = arfit(amp_ts, p_amp, p_amp);
        [~, A_phs, C_phs] = arfit(phs_ts, p_phs, p_phs);

        % Reshape coefficient arrays
        K = freq_amount;
        A_amp_mat = reshape(A_amp, K, K, p_amp);
        A_phs_mat = reshape(A_phs, K, K, p_phs);

        % Compute VAR PSDs via transfer matrix
        Z = exp(-1j*2*pi*f*sampling_interval); % vector length nfft/2+1
        H_amp = zeros(length(f), K, K);
        H_phs = zeros(length(f), K, K);
        for fi = 1:length(f)
            % Amplitude
            M = eye(K);
            for jj = 1:p_amp
                M = M - A_amp_mat(:,:,jj) * Z(fi)^(-jj);
            end
            H_amp(fi,:,:) = M \ eye(K);
            % Phase
            M2 = eye(K);
            for jj = 1:p_phs
                M2 = M2 - A_phs_mat(:,:,jj) * Z(fi)^(-jj);
            end
            H_phs(fi,:,:) = M2 \ eye(K);
        end

        % Assume H_amp: [F×K×K],  C_amp: [K×K],  acc_var_amp.(sev): [F×K×K]
        %        H_phs: [F×K×K], C_phs: [K×K],   acc_var_phs.(sev): [F×K×K]
        
        % Preallocate once (outside MC loop)
        S_amp = zeros(length(f),K,K);
        S_phs = zeros(length(f),K,K);
        
        for fi = 1:length(f)
            Hf_amp = squeeze(H_amp(fi,:,:));   % [K×K]
            S_amp(fi,:,:) = Hf_amp * C_amp * Hf_amp'; 
            
            Hf_phs = squeeze(H_phs(fi,:,:));   % [K×K]
            S_phs(fi,:,:) = Hf_phs * C_phs * Hf_phs';
        end
        
        % Now accumulate, applying the 2*(1/fs) and real() in one go:
        acc_var_amp.(sev) = acc_var_amp.(sev) + 2*(1/fs)*abs(S_amp);
        acc_var_phs.(sev) = acc_var_phs.(sev) + 2*(1/fs)*abs(S_phs);
    end
end

%% Average Monte Carlo and store in psd_comparison

psd_comparison = struct('freq',f,'amplitude',[],'total_phs',[]);
for i = 1:numel(severities)
    sev = severities{i};
    psd_comparison.amplitude.(sev).periodogram  = acc_per_amp.(sev)/num_realizations;
    psd_comparison.amplitude.(sev).var_psd      = acc_var_amp.(sev)/num_realizations;
    psd_comparison.total_phs.(sev).periodogram  = acc_per_phs.(sev)/num_realizations;
    psd_comparison.total_phs.(sev).var_psd      = acc_var_phs.(sev)/num_realizations;
end

%% Triple-frequency amplitude PSD comparison (3×3)
cmap_emp = winter(numel(severities));
cmap_var = autumn(numel(severities));

% Preallocate handles for legend (two lines per severity)
h = gobjects(numel(severities)*2,1);

figure('Position',[50,50,1200,700]);
tl = tiledlayout(K, K, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for p = 1:K
    for q = 1:K
        ax = nexttile(tl);
        hold(ax,'on');

        for k = 1:numel(severities)
            sev  = severities{k};
            per  = psd_comparison.amplitude.(sev).periodogram(:,p,q);
            varp = psd_comparison.amplitude.(sev).var_psd(:,p,q);

            if p==1 && q==1
                % Capture handles only once for the shared legend
                h(2*k-1) = plot(ax, f,10*log10(per), '--', ...
                    'Color',cmap_emp(k,:), ...
                    'DisplayName',sprintf('%s – Periodogram',sev));
                h(2*k)   = plot(ax, f,10*log10(varp), '-', ...
                    'Color',cmap_var(k,:), ...
                    'LineWidth',1.5, ...
                    'DisplayName',sprintf('%s – VAR PSD',sev));
            else
                plot(ax, f,10*log10(per), '--', 'Color',cmap_emp(k,:), ...
                    'DisplayName',sprintf('%s – Periodogram',sev));
                plot(ax, f,10*log10(varp), '-', 'Color',cmap_var(k,:), ...
                    'LineWidth',1.5, ...
                    'DisplayName',sprintf('%s – VAR PSD',sev));
            end
        end

        hold(ax,'off');
        ax.XScale   = 'log';
        ax.XLim     = [1e-4*fs, 0.4*fs];
        ax.FontSize = font_size;
        title(ax, sprintf('%s – %s', frequency_bands{p}, frequency_bands{q}));

        if p == K, xlabel(ax,'Norm. freq. (× 1/T_I) [Hz]'); end
        if q == 1, ylabel(ax,'PSD [dB/Hz]');   end
        grid(ax,'on');
    end
end

% Place a single legend in the east tile of the layout
lgd = legend(h, 'FontSize',font_size);
lgd.Layout.Tile = 'east';

% Export to vector PDF
exportgraphics(gcf, fullfile(fig_dir,...
    'periodogram_vs_ar_psd_amp_triplefreq_cpssm.pdf'), ...
    'ContentType','vector');


%% Triple-frequency phase PSD comparison (3×3)
% (Exactly the same pattern, swapping amplitude→total_phs)

h = gobjects(numel(severities)*2,1);

figure('Position',[50,50,1200,700]);
tl = tiledlayout(K, K, ...
    'TileSpacing','compact', ...
    'Padding','compact');

for p = 1:K
    for q = 1:K
        ax = nexttile(tl);
        hold(ax,'on');

        for k = 1:numel(severities)
            sev  = severities{k};
            per  = psd_comparison.total_phs.(sev).periodogram(:,p,q);
            varp = psd_comparison.total_phs.(sev).var_psd(:,p,q);

            if p==1 && q==1
                h(2*k-1) = plot(ax, f,10*log10(per), '--', ...
                    'Color',cmap_emp(k,:), ...
                    'DisplayName',sprintf('%s – Periodogram',sev));
                h(2*k)   = plot(ax, f,10*log10(varp), '-', ...
                    'Color',cmap_var(k,:), ...
                    'LineWidth',1.5, ...
                    'DisplayName',sprintf('%s – VAR PSD',sev));
            else
                plot(ax, f,10*log10(per), '--', 'Color',cmap_emp(k,:), ...
                    'DisplayName',sprintf('%s – Periodogram',sev));
                plot(ax, f,10*log10(varp), '-', 'Color',cmap_var(k,:), ...
                    'LineWidth',1.5, ...
                    'DisplayName',sprintf('%s – VAR PSD',sev));
            end
        end

        hold(ax,'off');
        ax.XScale   = 'log';
        ax.XLim     = [1e-4*fs, 0.5*fs];
        ax.FontSize = font_size;
        title(ax, sprintf('%s – %s', frequency_bands{p}, frequency_bands{q}));

        if p == K, xlabel(ax,'Norm. freq. (× 1/T_I) [Hz]'); end
        if q == 1, ylabel(ax,'PSD [dB (rad^2/Hz)]');   end
        grid(ax,'on');
    end
end

lgd = legend(h, 'FontSize',font_size);
lgd.Layout.Tile = 'east';

exportgraphics(gcf, fullfile(fig_dir,...
    'periodogram_vs_ar_psd_phs_triplefreq_cpssm.pdf'), ...
    'ContentType','vector');


%% Export CSVs per severity
for i = 1:numel(severities)
    sev = severities{i};
    % Amplitude
    P = psd_comparison.amplitude.(sev).periodogram;
    V = psd_comparison.amplitude.(sev).var_psd;
    T_amp = table(f, squeeze(P(:,1,1)), squeeze(P(:,2,2)), squeeze(P(:,3,3)), ...
        squeeze(P(:,1,2)), squeeze(P(:,1,3)), squeeze(P(:,2,3)), ...
        squeeze(V(:,1,1)), squeeze(V(:,2,2)), squeeze(V(:,3,3)), ...
        squeeze(V(:,1,2)), squeeze(V(:,1,3)), squeeze(V(:,2,3)), ...
        'VariableNames', { ...
          'Freq_Hz', 'Auto_L1','Auto_L2','Auto_L5', ...
          'Cross_L1_L2','Cross_L1_L5','Cross_L2_L5', ...
          'Var_Auto_L1','Var_Auto_L2','Var_Auto_L5', ...
          'Var_Cross_L1_L2','Var_Cross_L1_L5','Var_Cross_L2_L5' } );
    writetable(T_amp, fullfile(csv_dir, ['psd_amp_' sev '.csv']));
    % Phase
    Pp = psd_comparison.total_phs.(sev).periodogram;
    Vp = psd_comparison.total_phs.(sev).var_psd;
    T_phs = table(f, squeeze(Pp(:,1,1)), squeeze(Pp(:,2,2)), squeeze(Pp(:,3,3)), ...
        squeeze(Pp(:,1,2)), squeeze(Pp(:,1,3)), squeeze(Pp(:,2,3)), ...
        squeeze(Vp(:,1,1)), squeeze(Vp(:,2,2)), squeeze(Vp(:,3,3)), ...
        squeeze(Vp(:,1,2)), squeeze(Vp(:,1,3)), squeeze(Vp(:,2,3)), ...
        'VariableNames', { ...
          'Freq_Hz', 'Auto_L1','Auto_L2','Auto_L5', ...
          'Cross_L1_L2','Cross_L1_L5','Cross_L2_L5', ...
          'Var_Auto_L1','Var_Auto_L2','Var_Auto_L5', ...
          'Var_Cross_L1_L2','Var_Cross_L1_L5','Var_Cross_L2_L5' } );
    writetable(T_phs, fullfile(csv_dir, ['psd_phs_' sev '.csv']));
end

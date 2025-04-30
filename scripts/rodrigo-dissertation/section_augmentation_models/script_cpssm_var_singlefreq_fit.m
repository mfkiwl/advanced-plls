% Docstring
%
% Script to identify the goodness of fit of the VAR model for 
% single-frequency scintillation ampltiude and phase time series 
% genenrated by the compact phase-screen-based scintillation model (CPSSM) 
% using the ARFIT algorithm.
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
mc_runs    = 100;
min_order  = 1;
max_order  = 30;
optimal_orders_amp   = zeros(mc_runs,numel(severities));
optimal_orders_full_phase = zeros(mc_runs,numel(severities));
optimal_orders_diffractive_phase = zeros(mc_runs,numel(severities));

seed = 1;
for mc_idx = 1:mc_runs
    for i = 1:numel(severities)
        severity = severities{i};
        rng(seed);
        [scint_ts, refractive_phase] = get_tppsm_data(cpssm_params.(severity), 'seed', seed);
        amp_ts    = abs(scint_ts);
        % NOTE: `get_corrected_phase` is a function from the
        % 'gnss_scintillation_simulator' git submodule.
        propagated_field_phase_ts_training  = get_corrected_phase(scint_ts);
        diffractive_phase_ts_training = wrapToPi(propagated_field_phase_ts_training - refractive_phase);

        [~, A_amp]   = arfit(amp_ts,   min_order, max_order);
        [~, A_propagated_field_phase] = arfit(propagated_field_phase_ts_training, min_order, max_order);
        [~, A_diffractive_phase] = arfit(diffractive_phase_ts_training, min_order, max_order);

        optimal_orders_amp(mc_idx,i)   = size(A_amp,2)/size(A_amp,1);
        optimal_orders_full_phase(mc_idx,i) = size(A_propagated_field_phase,2)/size(A_propagated_field_phase,1);
        optimal_orders_diffractive_phase(mc_idx,i) = size(A_diffractive_phase,2)/size(A_diffractive_phase,1);
        seed = seed + 1;
    end
end

orders       = min_order:max_order;
counts_amp   = zeros(numel(orders),numel(severities));
counts_propagated_field_phase = zeros(numel(orders),numel(severities));
counts_diffractive_phase = zeros(numel(orders),numel(severities));
for i = 1:numel(severities)
    counts_amp(:,i)   = histcounts(optimal_orders_amp(:,i),   [orders, orders(end)+1]);
    counts_propagated_field_phase(:,i) = histcounts(optimal_orders_full_phase(:,i), [orders, orders(end)+1]);
    counts_diffractive_phase(:,i) = histcounts(optimal_orders_diffractive_phase(:,i), [orders, orders(end)+1]);
end

figure('Position',[100,100,1800,400]);
subplot(1,3,1);
bar(orders,counts_amp,'grouped');
xlabel('VAR Model Order'); ylabel('Frequency');
title('Optimal VAR Orders (Propagated Scintillation field Amplitude)');
legend(severities,'Location','best'); grid on;
subplot(1,3,2);
bar(orders,counts_propagated_field_phase,'grouped');
xlabel('VAR Model Order'); ylabel('Frequency');
title('Optimal VAR Orders (Propagated scintillation field Phase---Refractive + Diffractive)');
legend(severities,'Location','best'); grid on;
subplot(1,3,3);
bar(orders,counts_diffractive_phase,'grouped');
xlabel('VAR Model Order'); ylabel('Frequency');
title('Optimal VAR Orders (Diffractive phase---Propagated scintillation field phase - Refractive phase)');
legend(severities,'Location','best'); grid on;

%% Obtain residuals for most frequent orders
[~, highest_freq_idx_amp]   = max(counts_amp,[],1);
[~, highest_freq_idx_propagated_field_phase] = max(counts_propagated_field_phase,[],1);
[~, highest_freq_idx_diff_phase] = max(counts_diffractive_phase,[],1);

residuals = struct('amplitude',[],'prop_field_phase',[], 'diff_phase', []);

seed = 10;
for i = 1:numel(severities)
    severity = severities{i};
    rng(seed);

    % Training the VAR model
    [scint_ts, refractive_phase] = get_tppsm_data(cpssm_params.(severity), 'seed', seed);
    amp_ts_training    = abs(scint_ts);
    propagated_field_phase_ts_training  = get_corrected_phase(scint_ts);
    diffractive_phase_ts_training = wrapToPi(propagated_field_phase_ts_training - refractive_phase);

    ord_amp = orders(highest_freq_idx_amp(i));
    ord_prop_field_phase = orders(highest_freq_idx_propagated_field_phase(i));
    ord_diff_phase = orders(highest_freq_idx_diff_phase(i));

    [w_amp, A_amp]   = arfit(amp_ts_training, ord_amp, ord_amp);
    [w_prop_field_phase, A_prop_field_phase] = arfit(propagated_field_phase_ts_training, ord_prop_field_phase, ord_prop_field_phase);
    [w_diff_phase, A_diff_phase] = arfit(diffractive_phase_ts_training, ord_diff_phase, ord_diff_phase);

    % Obtaining the residuals
    % NOTE: Another generation seed is used herein to evaluate the reisudues 
    % correlation with statistical significance.
    rng(seed+1);
    [scint_ts, refractive_phase] = get_tppsm_data(cpssm_params.(severity), 'seed', seed + 1);
    amp_ts_for_residue    = abs(scint_ts);
    propagated_field_phase_ts_for_residue  = get_corrected_phase(scint_ts);
    diffractive_phase_ts_for_residue = wrapToPi(propagated_field_phase_ts_for_residue - refractive_phase);

    [~, res_amp]   = arres(w_amp,   A_amp,   amp_ts_for_residue,   60);
    [~, res_prop_field_phase]   = arres(w_prop_field_phase,   A_prop_field_phase,   propagated_field_phase_ts_for_residue,   60);
    [~, res_diff_phase]   = arres(w_diff_phase,   A_diff_phase,   diffractive_phase_ts_for_residue,   60);

    residuals.amplitude.(severity) = [NaN(ord_amp,1);   res_amp];
    residuals.prop_field_phase.(severity) = [NaN(ord_prop_field_phase,1);   res_prop_field_phase];
    residuals.diff_phase.(severity) = [NaN(ord_diff_phase,1);   res_diff_phase];

    seed = seed + 1;
end

%% Compute one-sided ACFs of residuals

% Amount of lags on the ACF
lags        = 20;
stem_width  = 1.5;
markers     = struct('Weak','o','Moderate','s','Strong','^');
acfs        = struct('amplitude',[],'prop_field_phase',[], 'diff_phase', []);

for i = 1:numel(severities)
    severity = severities{i};
    
    % Removing the NaNs of the amplitude and phase residuals
    amp_res  = residuals.amplitude.(severity);
    amp_res  = amp_res(~isnan(amp_res));
    prop_field_phase_res  = residuals.prop_field_phase.(severity);
    prop_field_phase_res  = prop_field_phase_res(~isnan(prop_field_phase_res));
    diff_phase_res  = residuals.diff_phase.(severity);
    diff_phase_res  = diff_phase_res(~isnan(diff_phase_res));

    acf_amp          = xcorr(amp_res, amp_res, lags, 'normalized');
    acf_prop_field_phase          = xcorr(prop_field_phase_res, prop_field_phase_res, lags, 'normalized');
    acf_diff_phase          = xcorr(diff_phase_res, diff_phase_res, lags, 'normalized');

    acfs.amplitude.(severity) = acf_amp(lags+1:end);
    acfs.prop_field_phase.(severity) = acf_prop_field_phase(lags+1:end);
    acfs.diff_phase.(severity) = acf_diff_phase(lags+1:end);
end

time_lag = (0:lags) * sampling_interval;
figure('Position',[100,100,1200,400]);
subplot(1,3,1);
hold on;
for i = 1:numel(severities)
    severity = severities{i};
    stem(time_lag, acfs.amplitude.(severity), 'LineWidth', stem_width, ...
         'Marker', markers.(severity), 'DisplayName', severity);
end
hold off;
xlabel('Time Lag [s]'); ylabel('ACF Amplitude');
title('Amplitude Residuals ACF');
legend('Location','best'); grid on;

subplot(1,3,2);
hold on;
for i = 1:numel(severities)
    severity = severities{i};
    stem(time_lag, acfs.prop_field_phase.(severity), 'LineWidth', stem_width, ...
         'Marker', markers.(severity), 'DisplayName', severity);
end
hold off;
xlabel('Time Lag [s]'); ylabel('ACF Phase');
title('Prop. field phase Residuals ACF');
legend('Location','best'); grid on;

subplot(1,3,3);
hold on;
for i = 1:numel(severities)
    severity = severities{i};
    stem(time_lag, acfs.diff_phase.(severity), 'LineWidth', stem_width, ...
         'Marker', markers.(severity), 'DisplayName', severity);
end
hold off;
xlabel('Time Lag [s]'); ylabel('ACF Phase');
title('Diffractive phase Residuals ACF');
legend('Location','best'); grid on;

%% Periodogram vs. VAR PSD comparison

% Number of points on the frequency domain
nfft = 2^16;
% Sampling frequency in Hz
fs   = 1/sampling_interval;

% Struct for pre-allocating the frequency support, the periodograms and 
% the VAR model's.
psd_comparison = struct( ...
  'freq',     [], ...
  'amplitude', struct('periodogram',[],'var_psd',[]), ...
  'prop_field_phase',     struct('periodogram',[],'var_psd',[]), ...
  'diff_phase',     struct('periodogram',[],'var_psd',[]) ...
);

seed = 1;
for i = 1:numel(severities)
    sev = severities{i};
    rng(seed);

    %---- 1) Generate & center signals ------------------------------------
    [scint_ts, refractive_phase] = get_tppsm_data(cpssm_params.(sev), 'seed', seed);
    amp_ts   = abs(scint_ts);
    prop_field_phase = get_corrected_phase(scint_ts);
    diff_phase = wrapToPi(prop_field_phase - refractive_phase);

    % center for correlation
    amp_ctr = amp_ts   - mean(amp_ts);
    prop_field_phase_ctr = prop_field_phase - mean(prop_field_phase);
    diff_phase_ctr = diff_phase - mean(diff_phase);

    %---- 2) Get the autocovariance sequence ------------------------------
    % returns c_full with lags = -(N-1):(N-1)
    [c_amp_full, lags] = xcorr(amp_ctr,   'biased');
    [c_prop_field_phase_full, ~   ] = xcorr(prop_field_phase_ctr,   'biased');
    [c_diff_phase_full, ~   ] = xcorr(diff_phase_ctr,   'biased');

    %---- 3) FFT of autocovariance sequence → two-sided PSD --
    % Amount of samples of the cross-correlation function (2 * N - 1)
    Ncorr      = numel(c_amp_full);
    % Amount of samples remaing to complete the fft size.
    pad_amount = nfft - Ncorr;
    C_amp_pad  = [c_amp_full; zeros(pad_amount,1)];
    C_prop_field_phase_pad  = [c_prop_field_phase_full; zeros(pad_amount,1)];
    C_diff_phase_pad  = [c_diff_phase_full; zeros(pad_amount,1)];

    % Get the fft of the padded autocovariance sequences
    S_amp_full = fft(C_amp_pad);
    S_prop_field_phase_full = fft(C_prop_field_phase_pad);
    S_diff_phase_full = fft(C_diff_phase_pad);

    % one-sided indices: {1 ... nfft/2+1}
    half = 1:(nfft/2+1);
    
    % Remove the remaining imaginary parts of the estimated periodogram.
    % NOTE: Here we multiply the periodogram with the sampling interval to
    % normalize it such that we have its unit in [dB/Hz]
    periodogram_amp = real(abs(S_amp_full(half))) * sampling_interval;
    periodogram_prop_field_phase = real(abs(S_prop_field_phase_full(half))) * sampling_interval;
    periodogram_diff_phase = real(abs(S_diff_phase_full(half))) * sampling_interval;

    % store freq once
    if isempty(psd_comparison.freq)
        % Parse the frequency support in Hz
        psd_comparison.freq = (half-1)*(fs/nfft);
    end

    %---- 4) Fit AR(p) & compute VAR-PSD (unchanged) ----------------------
    p_amp = orders(highest_freq_idx_amp(i));
    p_prop_field_phase = orders(highest_freq_idx_propagated_field_phase(i));
    p_diff_phase = orders(highest_freq_idx_diff_phase(i));
    [~, A_amp, C_amp] = arfit(amp_ts,   p_amp, p_amp);
    [~, A_prop_field_phase, C_prop_field_phase] = arfit(prop_field_phase, p_prop_field_phase, p_prop_field_phase);
    [~, A_diff_phase, C_diff_phase] = arfit(diff_phase, p_diff_phase, p_diff_phase);

    a_amp = [1, -reshape(A_amp,1,[])];
    a_prop_field_phase = [1, -reshape(A_prop_field_phase,1,[])];
    a_diff_phase = [1, -reshape(A_diff_phase,1,[])];

    % evaluate H at each frequency bin
    z = exp(-1j*2*pi* psd_comparison.freq * sampling_interval);
    H_amp  = zeros(size(z));
    H_prop_field_phase = zeros(size(z));
    H_diff_phase = zeros(size(z));

    for k = 0:p_amp
      H_amp = H_amp + a_amp(k+1) * z.^(-k);
    end
    H_amp = 1 ./ H_amp;

    for k = 0:p_prop_field_phase
      H_prop_field_phase = H_prop_field_phase + a_prop_field_phase(k+1) * z.^(-k);
    end
    H_prop_field_phase = 1 ./ H_prop_field_phase;

    for k = 0:p_diff_phase
      H_diff_phase = H_diff_phase + a_diff_phase(k+1) * z.^(-k);
    end
    H_diff_phase = 1 ./ H_diff_phase;

    S_var_amp = (C_amp/fs)  * abs(H_amp).^2;
    S_var_prop_field_phase = (C_prop_field_phase/fs)  * abs(H_prop_field_phase).^2;
    S_var_diff_phase = (C_diff_phase/fs)  * abs(H_diff_phase).^2;

    %---- 5) Stash results ------------------------------------------------
    psd_comparison.amplitude.(sev).periodogram = periodogram_amp;
    psd_comparison.amplitude.(sev).var_psd            = S_var_amp;
    psd_comparison.prop_field_phase.(sev).periodogram     = periodogram_prop_field_phase;
    psd_comparison.prop_field_phase.(sev).var_psd         = S_var_prop_field_phase;
    psd_comparison.diff_phase.(sev).periodogram     = periodogram_diff_phase;
    psd_comparison.diff_phase.(sev).var_psd         = S_var_diff_phase;
    
    seed = seed + 1;
end

%---- Plot ----------------------------------------------------------------
cmap_p = winter(numel(severities));  % periodogram lines
cmap_v = autumn(numel(severities));  % VAR-PSD     lines

figure('Position',[100,100,600,1200]);
F = psd_comparison.freq;

% Amplitude
subplot(3,1,1); hold on;
h_amp = gobjects(2*numel(severities),1);
for k=1:numel(severities)
  sev = severities{k};
  idx = 2*(k-1)+1;
  h_amp(idx)   = plot(F, 10*log10(psd_comparison.amplitude.(sev).periodogram), '--', ...
                      'Color',cmap_p(k,:), 'LineWidth',1, 'DisplayName',[sev ' – Corr. PSD']);
  h_amp(idx+1) = plot(F, 10*log10(psd_comparison.amplitude.(sev).var_psd),           '-', ...
                      'Color',cmap_v(k,:), 'LineWidth',2, 'DisplayName',[sev ' – VAR PSD']);
end
hold off;
set(gca,'XScale','log','XLim',[1e-4*fs,0.4*fs]);
xlabel('Frequency [Hz]'); ylabel('Power [dB]');
title('Amplitude: Corr-based PSD vs. VAR PSD');
grid on; legend(h_amp,'Location','best');

% Propagated field phase
subplot(3,1,2); hold on;
h_phs = gobjects(2*numel(severities),1);
for k=1:numel(severities)
  sev = severities{k};
  idx = 2*(k-1)+1;
  h_phs(idx)   = plot(F, 10*log10(psd_comparison.prop_field_phase.(sev).periodogram), '--', ...
                      'Color',cmap_p(k,:), 'LineWidth',1, 'DisplayName',[sev ' – Corr. PSD']);
  h_phs(idx+1) = plot(F, 10*log10(psd_comparison.prop_field_phase.(sev).var_psd),           '-', ...
                      'Color',cmap_v(k,:), 'LineWidth',2, 'DisplayName',[sev ' – VAR PSD']);
end
hold off;
set(gca,'XScale','log','XLim',[1e-4*fs,0.4*fs]);
xlabel('Frequency [Hz]'); ylabel('Power [dB]');
title('Phase: Corr-based PSD vs. VAR PSD');
grid on; legend(h_phs,'Location','best');

% Diffractive phase
subplot(3,1,3); hold on;
h_phs = gobjects(2*numel(severities),1);
for k=1:numel(severities)
  sev = severities{k};
  idx = 2*(k-1)+1;
  h_phs(idx)   = plot(F, 10*log10(psd_comparison.diff_phase.(sev).periodogram), '--', ...
                      'Color',cmap_p(k,:), 'LineWidth',1, 'DisplayName',[sev ' – Corr. PSD']);
  h_phs(idx+1) = plot(F, 10*log10(psd_comparison.diff_phase.(sev).var_psd),           '-', ...
                      'Color',cmap_v(k,:), 'LineWidth',2, 'DisplayName',[sev ' – VAR PSD']);
end
hold off;
set(gca,'XScale','log','XLim',[1e-4*fs,0.4*fs]);
xlabel('Frequency [Hz]'); ylabel('Power [dB]');
title('Phase: Corr-based PSD vs. VAR PSD');
grid on; legend(h_phs,'Location','best');
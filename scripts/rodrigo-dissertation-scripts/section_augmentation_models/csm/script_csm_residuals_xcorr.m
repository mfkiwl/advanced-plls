% analyze_csm_residual_acf.m

clearvars; clc;
addpath(genpath(fullfile(pwd,'..','..','..','..','libs')));

%% 1) PARAMETERS

mc_runs     = 2;        % Monte Carlo runs
lags        = 0:25;        % lags for ACF
nLags       = numel(lags);

% (a) your known “top 3” orders for each severity:
%    rows = [weak; moderate; strong], columns = [order1, order2, order3]
topOrders = [4, 5, 6;     % example for amplitude
             3, 4, 5];    % example for phase
% — replace the above with your actual three orders

severities = {'Weak','Moderate','Strong'};
csm_params.Weak     = struct('S4',0.2,'tau0',1,  'simulation_time',25,'sampling_interval',0.01);
csm_params.Moderate = struct('S4',0.5,'tau0',0.6,'simulation_time',25,'sampling_interval',0.01);
csm_params.Strong   = struct('S4',0.9,'tau0',0.2,'simulation_time',25,'sampling_interval',0.01);

%% 2) PREALLOCATE

% ACF storage: dims = [nLags × mc_runs × (3 orders) × (3 severities)]
acf_amp   = zeros(nLags, mc_runs, 3, 3);
acf_phase = zeros(nLags, mc_runs, 3, 3);

%% 3) MONTE CARLO LOOP

seed = 1;
for run = 1:mc_runs
    rng(seed);
    for si = 1:3                      % for each severity
        sev = severities{si};
        % generate one draw:
        ts = get_csm_data(csm_params.(sev));
        amp = abs(ts);
        ph  = atan2(imag(ts), real(ts));

        for kOrd = 1:3                 % for each of the top 3 orders
            % amplitude
            p = topOrders(1, kOrd);
            [wA,A_amp,~] = arfit( amp, p, p );
            resA = arres(wA, A_amp, amp, max(lags));
              % full ACF (lags –25:25)
            full_acf = xcorr( reshape(resA,[],1), max(lags), 'coeff' );
            center   = max(lags) + 1;             % index corresponding to lag 0
            % extract only lags 0:25
            acf_amp(:,run,kOrd,si) = full_acf(center + lags);
            % phase
            p = topOrders(2, kOrd);
            [wP,A_ph,~] = arfit( ph, p, p );
            resP = arres(wP, A_ph, ph, max(lags));
            full_acf = xcorr( reshape(resP,[],1), max(lags), 'coeff' );
            acf_phase(:,run,kOrd,si) = full_acf(center + lags);
        end
    end
    seed = seed + 1;
end

%% 4) COMPUTE MEAN & VAR

mean_amp   = squeeze( mean(acf_amp,   2) );   % [nLags×3×3]
var_amp    = squeeze( var (acf_amp,   0, 2) );

mean_phase = squeeze( mean(acf_phase, 2) );
var_phase  = squeeze( var (acf_phase, 0, 2) );

%% 5) PLOTTING

figure('Position',[100 100 1000 400]);

% --- Amplitude subplot ---
subplot(1,2,1); hold on;
colors = lines(9);
ci = 1;
for si = 1:3
  for kOrd = 1:3
    lbl = sprintf('%s–p=%d', severities{si}, topOrders(1,kOrd));
    plot(lags, mean_amp(:,kOrd,si), 'Color',colors(ci,:),'LineWidth',1.5);
    % shaded var ribbon
    yy = sqrt(var_amp(:,kOrd,si));
    fill([lags,fliplr(lags)], ...
         [mean_amp(:,kOrd,si)'+yy'; mean_amp(:,kOrd,si)'-yy'], ...
         colors(ci,:), 'FaceAlpha',0.1,'EdgeColor','none');
    ci = ci+1;
  end
end
title('Residual ACF — Amplitude');
xlabel('Lag'); ylabel('ACF');
legend(); grid on;  % fill in legend entries if you wish

% --- Phase subplot ---
subplot(1,2,2); hold on;
ci = 1;
for si = 1:3
  for kOrd = 1:3
    lbl = sprintf('%s–p=%d', severities{si}, topOrders(2,kOrd));
    plot(lags, mean_phase(:,kOrd,si), 'Color',colors(ci,:),'LineWidth',1.5);
    yy = sqrt(var_phase(:,kOrd,si));
    fill([lags,fliplr(lags)], ...
         [mean_phase(:,kOrd,si)'+yy'; mean_phase(:,kOrd,si)'-yy'], ...
         colors(ci,:), 'FaceAlpha',0.1,'EdgeColor','none');
    ci = ci+1;
  end
end
title('Residual ACF — Phase');
xlabel('Lag'); ylabel('ACF');
legend(); grid on;


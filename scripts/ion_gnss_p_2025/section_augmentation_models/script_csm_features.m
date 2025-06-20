%% script_csm_residuals_tests.m
clearvars; clc;
addpath(genpath(fullfile(pwd,'..','..','..','libs')));

%% Simulation parameters & severities
simulation_time   = 300;
sampling_interval = 0.01;
severities        = {'Weak','Moderate','Strong'};
csm_params = struct( ...
  'Weak',    struct('S4',0.2, 'tau0',1.0, 'simulation_time',simulation_time,'sampling_interval',sampling_interval),...
  'Moderate',struct('S4',0.5, 'tau0',0.6, 'simulation_time',simulation_time,'sampling_interval',sampling_interval),...
  'Strong',  struct('S4',0.9, 'tau0',0.2, 'simulation_time',simulation_time,'sampling_interval',sampling_interval) ...
);

%% Preallocate residuals struct
residues = struct('amplitude',struct(),'phase',struct());

%% First: compute "most frequent" AR orders (or hard-code them if known)
% Here we simply pick order = 8 for both amp & phase as an example; 
% replace with your logic for `highest_freq_idx_*` if you have it.
ord_amp   = 8;
ord_phase = 8;

for i = 1:numel(severities)
  sev = severities{i};
  rng(i);  % reproducibility
  data     = get_csm_data(csm_params.(sev));
  amp_ts   = abs(data);
  phs_ts   = atan2(imag(data),real(data));

  % fit AR and extract residuals
  [wA,Aamp] = arfit(amp_ts,ord_amp,ord_amp);
  [~,   rA] = arres(wA, Aamp, amp_ts, 20);
  [wP, Apha] = arfit(phs_ts,ord_phase,ord_phase);
  [~,   rP] = arres(wP, Apha, phs_ts, 20);

  % pad with NaNs to keep original length
  residues.amplitude.(sev) = [NaN(ord_amp,1); rA];
  residues.phase.(sev)     = [NaN(ord_phase,1); rP];
end

%% Define tests to run
tests   = {'ADF','KPSS','LMC','PP','VR','ARCH','LBQ'};
nTests  = numel(tests);
nSev    = numel(severities);

% storage: p-values, stats, critical values, VR‐ratios (strings)
p_amp    = nan(nTests,nSev);   stat_amp = nan(nTests,nSev);
crit_amp = nan(nTests,nSev);   vr_amp   = strings(nTests,nSev);
p_phs    = nan(nTests,nSev);   stat_phs = nan(nTests,nSev);
crit_phs = nan(nTests,nSev);   vr_phs   = strings(nTests,nSev);

%% A) Stationarity / Unit-Root Tests
%    – ADF, KPSS, Leybourne–McCabe, Phillips–Perron, Variance-Ratio
for s = 1:nSev
  % pick the residual series (drop NaNs)
  xA = residues.amplitude.(severities{s});
  xA = xA(~isnan(xA));
  xP = residues.phase.(severities{s});
  xP = xP(~isnan(xP));

  for t = 1:nTests
    name = tests{t};
    switch name
      case 'ADF'
        [h,p,stat,cVal] = adftest(xA,'Model','TS','Lags',0);
        vrA = '--';
        [~,p_p,stat_p,cVal_p] = adftest(xP,'Model','TS','Lags',0);
        vrP = '--';
      case 'KPSS'
        [h,p,stat,cVal] = kpsstest(xA,'Trend',true,'Lags',0);
        vrA = '--';
        [~,p_p,stat_p,cVal_p] = kpsstest(xP,'Trend',true,'Lags',0);
        vrP = '--';
      case 'LMC'
        [h,p,stat,cVal] = lmctest(xA,'Trend',true,'Lags',0);
        vrA = '--';
        [~,p_p,stat_p,cVal_p] = lmctest(xP,'Trend',true,'Lags',0);
        vrP = '--';
      case 'PP'
        [h,p,stat,cVal] = pptest(xA,'Lags',0,'Model','AR');
        vrA = '--';
        [~,p_p,stat_p,cVal_p] = pptest(xP,'Lags',0,'Model','AR');
        vrP = '--';
      case 'VR'
        [h,p,stat,cVal,ratio] = vratiotest(xA,'Period',2,'IID',false);
        vrA    = ratio;
        [~,p_p,stat_p,cVal_p,ratio_p] = vratiotest(xP,'Period',2,'IID',false);
        vrP    = ratio_p;
      case 'ARCH'
        [h,p,stat,cVal]   = archtest(xA,'Lags',1);
        vrA    = '--';
        [~,p_p,stat_p,cVal_p] = archtest(xP,'Lags',1);
        vrP    = '--';
      case 'LBQ'
        [h,p,stat,cVal]   = lbqtest(xA,'Lags',20);
        vrA    = '--';
        [~,p_p,stat_p,cVal_p] = lbqtest(xP,'Lags',20);
        vrP    = '--';
    end

    % store amplitude
    p_amp(t,s)    = p;
    stat_amp(t,s) = stat;
    crit_amp(t,s) = cVal;
    vr_amp(t,s)   = vrA;

    % store phase
    p_phs(t,s)    = p_p;
    stat_phs(t,s) = stat_p;
    crit_phs(t,s) = cVal_p;
    vr_phs(t,s)   = vrP;
  end
end

%% Build & display result tables

% Column names: Test, then for each sev → {p,stat,crit,vr}
varNames = ['Test', reshape([strcat(severities,'_p'); strcat(severities,'_stat'); ...
                             strcat(severities,'_crit'); strcat(severities,'_vr')],1,[])];

% Fill cell array
AmpCell = cell(nTests,1+4*nSev);
PhsCell = cell(nTests,1+4*nSev);
for t=1:nTests
  AmpCell{t,1} = tests{t};
  PhsCell{t,1} = tests{t};
  for s=1:nSev
    col0 = 1 + (s-1)*4;
    AmpCell{t,col0+1} = sprintf('%.4f',   p_amp(t,s));
    AmpCell{t,col0+2} = sprintf('%.4f',   stat_amp(t,s));
    AmpCell{t,col0+3} = sprintf('%.4f',   crit_amp(t,s));
    AmpCell{t,col0+4} = vr_amp(t,s);
    PhsCell{t,col0+1} = sprintf('%.4f',   p_phs(t,s));
    PhsCell{t,col0+2} = sprintf('%.4f',   stat_phs(t,s));
    PhsCell{t,col0+3} = sprintf('%.4f',   crit_phs(t,s));
    PhsCell{t,col0+4} = vr_phs(t,s);
  end
end

T_amp = cell2table(AmpCell,'VariableNames',varNames);
T_phs = cell2table(PhsCell,'VariableNames',varNames);

fprintf('\n=== A) Stationarity / Unit-Root Tests [and B & C mixed in] ===\n');
disp('--- Amplitude residuals ---');
disp(T_amp);
disp('--- Phase residuals ---');
disp(T_phs);

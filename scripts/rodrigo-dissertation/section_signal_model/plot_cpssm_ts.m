clearvars; clc;
addpath(genpath(fullfile(pwd,"..","..","..")));

% make results folder
if ~exist("results_pdf","dir")
    mkdir("results_pdf");
end

%% parameters
sim_time        = 150;
t_samp          = 0.01;
time            = (t_samp:t_samp:sim_time).';
zoom_idx        = 10000:15000;
time_zoom       = time(zoom_idx);

seed            = 10; rng(seed);

severities      = {'Weak','Moderate','Strong'};
frequency_bands = {'L1','L2','L5'};
freq_amount     = numel(frequency_bands);

% CPSSM calls
cpssm_params = struct(...
  'weak',     {'weak',     'is_enable_cmd_print',false,'simulation_time',sim_time,'sampling_interval',t_samp,'rhof_veff_ratio',1.5},...
  'moderate', {'moderate', 'is_enable_cmd_print',false,'simulation_time',sim_time,'sampling_interval',t_samp,'rhof_veff_ratio',0.8},...
  'strong',   {'strong',   'is_enable_cmd_print',false,'simulation_time',sim_time,'sampling_interval',t_samp,'rhof_veff_ratio',0.27}...
);

[weak_scint_ts,      weak_refr_phase_ts]     = get_tppsm_multifreq_data(cpssm_params.weak,     'seed',seed);
[moderate_scint_ts,  moderate_refr_phase_ts] = get_tppsm_multifreq_data(cpssm_params.moderate, 'seed',seed);
[strong_scint_ts,    strong_refr_phase_ts]   = get_tppsm_multifreq_data(cpssm_params.strong,   'seed',seed);

% intensities (dB)
weak_int_dB   = 10*log10(abs(weak_scint_ts    ).^2);
mod_int_dB    = 10*log10(abs(moderate_scint_ts).^2);
strong_int_dB = 10*log10(abs(strong_scint_ts  ).^2);

% total phase via get_corrected_phase
weak_tot   = [ get_corrected_phase(weak_scint_ts(:,1)),     get_corrected_phase(weak_scint_ts(:,2)),     get_corrected_phase(weak_scint_ts(:,3)) ];
mod_tot    = [ get_corrected_phase(moderate_scint_ts(:,1)), get_corrected_phase(moderate_scint_ts(:,2)), get_corrected_phase(moderate_scint_ts(:,3)) ];
strong_tot = [ get_corrected_phase(strong_scint_ts(:,1)),   get_corrected_phase(strong_scint_ts(:,2)),   get_corrected_phase(strong_scint_ts(:,3)) ];

% diffractive = total − refractive
weak_diff   = weak_tot   - weak_refr_phase_ts;
mod_diff    = mod_tot    - moderate_refr_phase_ts;
strong_diff = strong_tot - strong_refr_phase_ts;

% wrapped diffractive
weak_diff_w   = wrapToPi(weak_diff);
mod_diff_w    = wrapToPi(mod_diff);
strong_diff_w = wrapToPi(strong_diff);

% bundle data by severity for easy indexing
intensity_dB = { weak_int_dB,   mod_int_dB,   strong_int_dB   };
total_phase  = { weak_tot,      mod_tot,      strong_tot      };
refr_phase   = { weak_refr_phase_ts, moderate_refr_phase_ts, strong_refr_phase_ts };
diff_phase   = { weak_diff,     mod_diff,     strong_diff     };
diff_phase_w = { weak_diff_w,   mod_diff_w,   strong_diff_w   };

% plotting styles
cmap               = lines(3);            % one color per band
legFreq            = frequency_bands;     % {'L1','L2','L5'}
plot_order         = [3 2 1];             % L5, L2, L1
font_size          = 13;
legend_font_size   = 8.5;
line_width         = 1.5;

%% loop over each severity
for s = 1:3
  sev = severities{s};

  %% Full‐time figure
  fig = figure('Position',[0,0,1000,600]);
  tiledlayout(4,1,'TileSpacing','tight');

  % 1) Intensity [dB]
  ax = nexttile;
  hold(ax,'on');
  for idx = plot_order
    plot(ax, time, intensity_dB{s}(:,idx), 'Color',cmap(idx,:), 'LineWidth',line_width);
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$10\log_{10}(\rho_{l}^2[k])$ [dB]','Interpreter','latex');
  legend(ax, legFreq(plot_order), 'Location','southwest', 'FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % 2) Total (­) + Refractive (:) Phase [rad]
  ax = nexttile;
  hold(ax,'on');
  labels = {};
  for idx = plot_order
    plot(ax, time, total_phase{s}(:,idx), 'Color',cmap(idx,:), 'LineWidth',line_width,'LineStyle','-');
    plot(ax, time, refr_phase{s}(:,idx),  'Color',cmap(idx,:), 'LineWidth',line_width,'LineStyle',':');
    labels{end+1} = sprintf('%s total',     legFreq{idx});
    labels{end+1} = sprintf('%s refractive',legFreq{idx});
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$\phi_{\mathrm{I},l}[k]$ [rad]','Interpreter','latex');
  legend(ax, labels, 'Location','southwest', 'FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % 3) Diffractive (unwrapped) [rad]
  ax = nexttile;
  hold(ax,'on');
  for idx = plot_order
    plot(ax, time, diff_phase{s}(:,idx), 'Color',cmap(idx,:), 'LineWidth',line_width);
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$\phi_{\mathrm{D},l}[k]$ [rad]','Interpreter','latex');
  legend(ax, legFreq(plot_order), 'Location','southwest','FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % 4) Diffractive (wrapped) [rad]
  ax = nexttile;
  hold(ax,'on');
  for idx = plot_order
    plot(ax, time, diff_phase_w{s}(:,idx), 'Color',cmap(idx,:), 'LineWidth',line_width);
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$\phi_{\mathrm{D},l}[k]$ [rad]','Interpreter','latex');
  xlabel(ax,'Time [s]','Interpreter','latex');
  legend(ax, legFreq(plot_order), 'Location','southwest','FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % save full‐time PDF
  exportgraphics(fig, fullfile('results_pdf',sprintf('cpssm_%s_full.pdf',lower(sev))), 'ContentType','vector');


  %% Zoomed figure
  fig = figure('Position',[0,0,1000,600]);
  tiledlayout(4,1,'TileSpacing','tight');

  % 1) Intensity zoom
  ax = nexttile;
  hold(ax,'on');
  for idx = plot_order
    plot(ax, time_zoom, intensity_dB{s}(zoom_idx,idx), 'Color',cmap(idx,:), 'LineWidth',line_width);
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$10\log_{10}(\rho_{l}^2[k])$ [dB]','Interpreter','latex');
  legend(ax, legFreq(plot_order), 'Location','southwest','FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % 2) Total + Refractive zoom
  ax = nexttile;
  hold(ax,'on');
  labels = {};
  for idx = plot_order
    plot(ax, time_zoom, total_phase{s}(zoom_idx,idx), 'Color',cmap(idx,:), 'LineWidth',line_width,'LineStyle','-');
    plot(ax, time_zoom, refr_phase{s}(zoom_idx,idx),  'Color',cmap(idx,:), 'LineWidth',line_width,'LineStyle',':');
    labels{end+1} = sprintf('%s total',     legFreq{idx});
    labels{end+1} = sprintf('%s refractive',legFreq{idx});
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$\phi_{\mathrm{I},l}[k]$ [rad]','Interpreter','latex');
  legend(ax, labels, 'Location','southwest','FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % 3) Diffractive unwrapped zoom
  ax = nexttile;
  hold(ax,'on');
  for idx = plot_order
    plot(ax, time_zoom, diff_phase{s}(zoom_idx,idx), 'Color',cmap(idx,:), 'LineWidth',line_width);
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$\phi_{\mathrm{D},l}[k]$ [rad]','Interpreter','latex');
  legend(ax, legFreq(plot_order), 'Location','southwest','FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % 4) Diffractive wrapped zoom
  ax = nexttile;
  hold(ax,'on');
  for idx = plot_order
    plot(ax, time_zoom, diff_phase_w{s}(zoom_idx,idx), 'Color',cmap(idx,:), 'LineWidth',line_width);
  end
  hold(ax,'off');
  grid(ax,'on'); grid(ax,'minor');
  ylabel(ax,'$\phi_{\mathrm{D},l}[k]$ [rad]','Interpreter','latex');
  xlabel(ax,'Time [s]','Interpreter','latex');
  legend(ax, legFreq(plot_order), 'Location','southwest','FontSize',legend_font_size, 'Direction','reverse');
  set(ax,'FontSize',font_size,'FontName','Times New Roman');

  % save zoomed PDF
  exportgraphics(fig, fullfile('results_pdf',sprintf('cpssm_%s_zoom.pdf',lower(sev))), 'ContentType','vector');

end

clearvars; clc;
addpath(genpath(fullfile(pwd,"..","..","..")));

% make results folder
if ~exist("results","dir")
    mkdir("results");
end

%% parameters
sim_time        = 150;
t_samp          = 0.01;
time            = (t_samp:t_samp:sim_time).';
zoom_idx        = 10000:15000;
time_zoom       = time(zoom_idx);

seed            = 10; rng(seed);

severities      = {'Weak','Moderate','Strong'};
frequency_bands = {'L1'};
freq_amount     = numel(frequency_bands);

% CPSSM calls
cpssm_params = struct(...
  'weak',     {'weak',     'is_enable_cmd_print',false,'simulation_time',sim_time,'sampling_interval',t_samp,'rhof_veff_ratio',1.5},...
  'strong',   {'strong',   'is_enable_cmd_print',false,'simulation_time',sim_time,'sampling_interval',t_samp,'rhof_veff_ratio',0.27}...
);

[weak_scint_ts, weak_phi_R] = get_tppsm_multifreq_data(cpssm_params.weak, 'seed',seed);
[strong_scint_ts, strong_phi_R] = get_tppsm_multifreq_data(cpssm_params.strong, 'seed',seed);

% intensities (dB)
weak_int_dB   = 10*log10(abs(weak_scint_ts(:,1)).^2);
strong_int_dB = 10*log10(abs(strong_scint_ts(:,1)).^2);

% total phase via get_corrected_phase
weak_phi_I = get_corrected_phase(weak_scint_ts(:,1));
strong_phi_I = get_corrected_phase(strong_scint_ts(:,1));

% diffractive = total − refractive
weak_diff   = weak_phi_I - weak_phi_R(:,1);
strong_diff = strong_phi_I - strong_phi_R(:,1);

% wrapped diffractive
weak_diff_w   = wrapToPi(weak_diff);
strong_diff_w = wrapToPi(strong_diff);

% plotting styles
cmap               = lines(2); % one color per band
legFreq            = frequency_bands; % {'L1','L2','L5'}
plot_order         = [3 2 1]; % L5, L2, L1
font_size          = 13;
legend_font_size   = 8.5;
line_width         = 1.5;

%% Generate plots
figure('Position',[50,50,1000,600]);
tiledlayout(3,1,"TileSpacing","tight");

font_size = 13;
line_width = 1.5;
nexttile;
hold on;
plot(time,strong_int_dB, 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,weak_int_dB, 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$10 \mathrm{log}10(\rho^2 [k])$ [dB]','interpreter', 'latex');
legend({'Strong', 'Weak'}, 'location', 'southwest');
set(gca, 'FontSize', font_size, 'fontname', 'Times New Roman');

nexttile;
hold on;
plot(time,strong_diff_w, 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,weak_diff_w, 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\phi_{\mathrm{D}}[k]$ [rad]','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
set(gca, 'FontSize', font_size, 'fontname', 'Times New Roman');
exportgraphics(gcf, 'results/csm_ts.pdf', 'ContentType', 'vector');

nexttile;
hold on;
plot(time,strong_diff, 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,weak_diff, 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\mathrm{UW}\{\phi_{\mathrm{D}}[k]\}$ [rad]','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
set(gca, 'FontSize', font_size, 'fontname', 'Times New Roman');
exportgraphics(gcf, 'results/cpssm_ts.pdf', 'ContentType', 'vector');

% nexttile;
% hold on;
% plot(time,strong_phi_I, 'LineWidth',line_width, 'Color', cmap(2,:));
% plot(time,weak_phi_I, 'LineWidth',line_width, 'Color', cmap(1,:));
% hold off;
% grid on;
% grid("minor");
% ylabel('$\phi_{\mathrm{I}} [k]$ [rad]','Interpreter','latex');
% set(gca, 'FontSize', font_size, 'fontname', 'Times New Roman');

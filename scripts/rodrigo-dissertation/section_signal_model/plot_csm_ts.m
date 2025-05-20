clearvars; clc;

addpath(genpath(fullfile(pwd, "..", "..", "..")));

if ~isempty(fullfile(pwd,'results_pdf'))
    mkdir('results_pdf');
end

sim_time = 150;
t_samp = 0.01;

seed = 2;
rng(seed);

csm_params_weak     = struct('S4', 0.2, 'tau0', 1.0, 'simulation_time', sim_time, 'sampling_interval', t_samp);
csm_params_moderate = struct('S4', 0.5, 'tau0', 0.6, 'simulation_time', sim_time, 'sampling_interval', t_samp);
csm_params_strong   = struct('S4', 0.9, 'tau0', 0.2, 'simulation_time', sim_time, 'sampling_interval', t_samp);

weak_ts = get_csm_data(csm_params_weak);
moderate_ts = get_csm_data(csm_params_moderate);
strong_ts = get_csm_data(csm_params_strong);

time = t_samp:t_samp:sim_time;
zoom_idx = 10000:15000;
time_zoomed = time(zoom_idx);

cmap = lines(3);

%% Generate plots
figure('Position',[50,50,1000,600]);
tiledlayout(3,1,"TileSpacing","tight");

font_size = 13;
line_width = 1.5;
nexttile;
hold on;
plot(time,10*log10(abs(strong_ts).^2), 'LineWidth',line_width, 'Color', cmap(3,:));
plot(time,10*log10(abs(moderate_ts).^2), 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,10*log10(abs(weak_ts).^2), 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$10 \mathrm{log}10(\rho_1^2 [k])$','interpreter', 'latex');
legend({'Strong', 'Moderate', 'Weak'}, 'location', 'southwest');
set(gca, 'FontSize', font_size);

nexttile;
hold on;
plot(time,angle(strong_ts), 'LineWidth',line_width, 'Color', cmap(3,:));
plot(time,angle(moderate_ts), 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,angle(weak_ts), 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\phi_{\mathrm{I},1}[k]$ [rad]','Interpreter','latex');
set(gca, 'FontSize', font_size);

nexttile;
hold on;
plot(time,phase(strong_ts), 'LineWidth',line_width, 'Color', cmap(3,:));
plot(time,phase(moderate_ts), 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,phase(weak_ts), 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\mathrm{unwrap}\{\phi_{I,1}[k]\}$ [rad]','Interpreter','latex');
xlabel('Time [s]');
set(gca, 'FontSize', font_size);
exportgraphics(gcf, 'results_pdf/csm_ts.pdf');

%% Generate plots - zoomed
figure('Position',[50,50,1000,600]);
tiledlayout(3,1,"TileSpacing","tight");

font_size = 14;
line_width = 1.5;
nexttile;
hold on;
plot(time_zoomed,10*log10(abs(strong_ts(zoom_idx)).^2), 'LineWidth',line_width, 'Color', cmap(3,:));
plot(time_zoomed,10*log10(abs(moderate_ts(zoom_idx)).^2), 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time_zoomed,10*log10(abs(weak_ts(zoom_idx)).^2), 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$10 \mathrm{log}10(\rho_1^2 [k])$','interpreter', 'latex');
legend({'Strong', 'Moderate', 'Weak'}, 'location', 'southwest');
set(gca, 'FontSize', font_size);

nexttile;
hold on;
plot(time_zoomed,angle(strong_ts(zoom_idx)), 'LineWidth',line_width, 'Color', cmap(3,:));
plot(time_zoomed,angle(moderate_ts(zoom_idx)), 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time_zoomed,angle(weak_ts(zoom_idx)), 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\phi_{\mathrm{I},1}[k]$ [rad]','Interpreter','latex');
set(gca, 'FontSize', font_size);

nexttile;
hold on;
plot(time_zoomed,phase(strong_ts(zoom_idx)), 'LineWidth',line_width, 'Color', cmap(3,:));
plot(time_zoomed,phase(moderate_ts(zoom_idx)), 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time_zoomed,phase(weak_ts(zoom_idx)), 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\mathrm{unwrap}\{\phi_{I,1}[k]\}$ [rad]','Interpreter','latex');
xlabel('Time [s]');
set(gca, 'FontSize', font_size);
exportgraphics(gcf, 'results_pdf/csm_ts_zoom.pdf');

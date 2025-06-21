clearvars; clc;

addpath(genpath(fullfile(pwd, "..", "..", "..")));

if ~exist("results","dir")
    mkdir('results');
end

sim_time = 150;
t_samp = 0.01;

seed = 5;
rng(seed);

csm_params_weak     = struct('S4', 0.2, 'tau0', 1.0, 'simulation_time', sim_time, 'sampling_interval', t_samp);
csm_params_strong   = struct('S4', 0.9, 'tau0', 0.2, 'simulation_time', sim_time, 'sampling_interval', t_samp);

weak_ts = get_csm_data(csm_params_weak);
strong_ts = get_csm_data(csm_params_strong);

weak_intensity_ts = 10*log10(abs(weak_ts).^2);
strong_intensity_ts = 10*log10(abs(strong_ts).^2);

weak_phase_ts = angle(weak_ts);
strong_phase_ts = angle(strong_ts);

weak_unwrapped_phase_ts = phase(weak_ts);
strong_unwrapped_phase_ts = phase(strong_ts);

time = t_samp:t_samp:sim_time;

cmap = lines(2);

%% Generate plots
figure('Position',[50,50,1000,600]);
tiledlayout(3,1,"TileSpacing","tight");

font_size = 13;
line_width = 1.5;
nexttile;
hold on;
plot(time,strong_intensity_ts, 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,weak_intensity_ts, 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$10 \mathrm{log}10(\rho^2 [k])$ [dB]','interpreter', 'latex');
legend({'Strong', 'Weak'}, 'location', 'southwest');
set(gca, 'FontSize', font_size, 'fontname', 'Times New Roman');

nexttile;
hold on;
plot(time,strong_phase_ts, 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,weak_phase_ts, 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\phi_{\mathrm{D}}[k]$ [rad]','Interpreter','latex');
set(gca, 'FontSize', font_size, 'fontname', 'Times New Roman');

nexttile;
hold on;
plot(time,strong_unwrapped_phase_ts, 'LineWidth',line_width, 'Color', cmap(2,:));
plot(time,weak_unwrapped_phase_ts, 'LineWidth',line_width, 'Color', cmap(1,:));
hold off;
grid on;
grid("minor");
ylabel('$\mathrm{UW}\{\phi_{\mathrm{D}}[k]\}$ [rad]','Interpreter','latex');
xlabel('Time [s]','Interpreter','latex');
set(gca, 'FontSize', font_size, 'fontname', 'Times New Roman');
exportgraphics(gcf, 'results/csm_ts.pdf', 'ContentType', 'vector');

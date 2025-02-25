function plot_KFAR_estimates_comparison(...
    time_vector, ...
    los_phase, ...
    rx_sig_csm, ...
    state_estimates_csm, ...
    state_estimates_none_under_csm, ...
    rx_sig_tppsm, ...
    state_estimates_tppsm, ...
    state_estimates_none_under_tppsm, ...
    doppler_profile, ...
    cnr, ...
    seed, ...
    process_noise_variance, ...
    ar_order)
% plot_csm_tppsm_time_series_comparison Plots a detailed time series comparison 
% between KF-AR and standard KF estimates under both CSM and TPPSM scenarios.
%
% The function creates a 2x1 tiled layout with the following plots:
%
%   Tile 1 (CSM | ):
%     - CSM + AWGN signal phase (detrended): unwrap(angle(rx_sig_csm))-los_phase
%     - KF-AR estimated LOS phase (detrended): state_estimates_csm(:,1)-los_phase
%     - KF-AR estimated scintillation phase: state_estimates_csm(:,length(doppler_profile)+1)
%     - Standard KF estimated LOS phase (detrended): state_estimates_none_under_csm(:,1)-los_phase
%
%   Tile 2 (TPPSM):
%     - TPPSM + AWGN signal phase (detrended): unwrap(angle(rx_sig_tppsm))-los_phase
%     - KF-AR estimated LOS phase (detrended): state_estimates_tppsm(:,1)-los_phase
%     - KF-AR estimated scintillation phase: state_estimates_tppsm(:,length(doppler_profile)+1)
%     - Standard KF estimated LOS phase (detrended): state_estimates_none_under_tppsm(:,1)-los_phase
%
% A global title is added along with an annotation showing the simulation parameters:
%   - Carrier-to-noise ratio (C/N0)
%   - Seed
%   - Process noise variance (for the discrete Wiener model)
%   - AR model order
%
% Finally, the figure is saved in both .fig and .pdf formats in a "figures" folder
% located one level above the current scripts folder.
%
% Inputs:
%   time_vector                   - [Nx1] vector with time stamps in seconds.
%   los_phase                     - [Nx1] vector with LOS phase (for detrending).
%   rx_sig_csm                  - [Nx1] complex received signal for the CSM scenario.
%   state_estimates_csm           - [NxM] matrix with KF-AR estimates for CSM.
%                                 First column is LOS phase; column (length(doppler_profile)+1)
%                                 is the AR (scintillation) phase estimate.
%   state_estimates_none_under_csm- [NxM] matrix with standard KF (none) estimates for CSM.
%   rx_sig_tppsm                  - [Nx1] complex received signal for the TPPSM scenario.
%   state_estimates_tppsm         - [NxM] matrix with KF-AR estimates for TPPSM.
%   state_estimates_none_under_tppsm - [NxM] matrix with standard KF (none) estimates for TPPSM.
%   doppler_profile               - vector; used to select the AR phase column as length(doppler_profile)+1.
%   cnr                           - Carrier-to-noise ratio (in dBHz).
%   seed                          - Seed used in the simulation.
%   process_noise_variance        - Process noise variance used in the discrete Wiener model.
%   ar_order                      - Order of the AR model (e.g., 3).
%
% Author: [Your Name]
% Date: [Today's Date]

% Create figure and tiled layout (2 rows x 1 column)
fig = figure('Position',[50,50,1200,800]);
tlo = tiledlayout(2,1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% Create a text string for simulation parameters
simParamsText = sprintf(['C/N0: %d dBHz; ' ...
                         'Seed: %d; ' ...
                         'Process Noise Variance: %g; ' ...
                         'AR Model Order: %d;'], cnr, seed, process_noise_variance, ar_order);
                     

%% --- Tile 1: CSM Scenario ---
nexttile;
% Compute phases (detrended by subtracting los_phase)
phase_csm          = unwrap(angle(rx_sig_csm)) - los_phase;
kf_ar_los_csm      = state_estimates_csm(:,1) - los_phase;
kf_ar_scint_csm    = state_estimates_csm(:, length(doppler_profile)+1);
kf_ar_joint_csm    = kf_ar_los_csm + kf_ar_scint_csm;
std_kf_csm         = state_estimates_none_under_csm(:,1) - los_phase;

% Plot each line with contrasting colors
hold on;
h1 = plot(time_vector, phase_csm, 'LineWidth', 1.5, 'Color', 'k');
h2 = plot(time_vector, kf_ar_los_csm, 'LineWidth', 1.5, 'Color', 'r');
h3 = plot(time_vector, kf_ar_scint_csm, 'LineWidth', 1.5, 'Color', 'g');
h4 = plot(time_vector, kf_ar_joint_csm, 'LineWidth', 1.5, 'Color', 'b');
h5 = plot(time_vector, std_kf_csm, 'LineWidth', 1.5, 'Color', 'magenta');
hold off;
grid on;
xlabel('Time [sec]');
ylabel('Phase [rad]');
title('CSM Severe Scenario');
legend({...
    'CSM + AWGN signal phase (detrended)',...
    'KF-AR estimated LOS phase (detrended)',...
    'KF-AR estimated scintillation phase',...
    'KF-AR estimated joint phase (LOS + scint, detrended)',...
    'Standard KF estimated LOS phase (detrended)'},...
    'Location', 'best');

%% --- Tile 2: TPPSM Scenario ---
nexttile;
% Compute phases for TPPSM (using same los_phase for detrending)
phase_tppsm          = unwrap(angle(rx_sig_tppsm)) - los_phase;
kf_ar_los_tppsm      = state_estimates_tppsm(:,1) - los_phase;
kf_ar_scint_tppsm    = state_estimates_tppsm(:, length(doppler_profile)+1);
kf_ar_joint_tppsm    = kf_ar_los_tppsm + kf_ar_scint_tppsm;
std_kf_tppsm       = state_estimates_none_under_tppsm(:,1) - los_phase;

hold on;
h1 = plot(time_vector, phase_tppsm, 'LineWidth', 1.5, 'Color', 'k');
h2 = plot(time_vector, kf_ar_los_tppsm, 'LineWidth', 1.5, 'Color', 'r');
h3 = plot(time_vector, kf_ar_scint_tppsm, 'LineWidth', 1.5, 'Color', 'g');
h4 = plot(time_vector, kf_ar_joint_tppsm, 'LineWidth', 1.5, 'Color', 'b');
h5 = plot(time_vector, std_kf_tppsm, 'LineWidth', 1.5, 'Color', 'magenta');
hold off;
grid on;
xlabel('Time [sec]');
ylabel('Phase [rad]');
title('TPPSM Severe Scenario');
legend({...
    'TPPSM + AWGN signal phase (detrended)',...
    'KF-AR estimated LOS phase (detrended)',...
    'KF-AR estimated scintillation phase',...
    'KF-AR estimated joint phase (LOS + scint, detrended)',...
    'Standard KF estimated LOS phase (detrended)'},...
    'Location', 'best');

% Global title and simulation parameters annotation
sgtitle(simParamsText);

%% Save the figure to disk
% Determine the figures folder relative to the current scripts folder.
% (Assuming current folder is ".../kalman_pll_testbench/scripts")
current_folder = pwd;
[baseFolder, ~, ~] = fileparts(current_folder);
figures_dir = fullfile(baseFolder, 'figures');
if ~exist(figures_dir, 'dir')
    mkdir(figures_dir);
end

% Adjust paper size for correct PDF export (in centimeters)
set(fig, 'PaperUnits', 'centimeters');
fig_pos = get(fig, 'Position'); % [left bottom width height] in pixels
fig_width_cm = fig_pos(3) * 2.54 / 96;   % conversion factor: pixels to cm (assuming 96 dpi)
fig_height_cm = fig_pos(4) * 2.54 / 96;
set(fig, 'PaperSize', [fig_width_cm fig_height_cm]);
set(fig, 'PaperPosition', [0 0 fig_width_cm fig_height_cm]);

% Save as high-quality PDF
fig_filename_pdf = fullfile(figures_dir, 'KF_AR_vs_StandardKF.pdf');
print(fig, fig_filename_pdf, '-dpdf', '-vector');

% Save as MATLAB figure (.fig)
fig_filename_fig = fullfile(figures_dir, 'KF_AR_vs_StandardKF.fig');
savefig(fig, fig_filename_fig);

end

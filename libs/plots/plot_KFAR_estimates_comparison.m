function plot_KFAR_estimates_comparison(...
    time_vector, ...
    los_phase, ...
    rx_sig_csm, ...
    kf_ar_csm, ...         % KF-AR estimates (no adaptive update) for CSM
    akf_ar_csm, ...        % AKF-AR estimates (simplified, no HL) for CSM
    ahl_kf_ar_csm, ...     % AHL-KF-AR estimates (simplified, HL) for CSM
    kf_std_csm, ...        % KF-std estimates (no adaptive update) for CSM
    akf_std_csm, ...       % AKF-std estimates (simplified, no HL) for CSM
    ahl_kf_std_csm, ...    % AHL-KF-std estimates (simplified, HL) for CSM
    rx_sig_tppsm, ...
    kf_ar_tppsm, ...       % KF-AR estimates for TPPSM
    akf_ar_tppsm, ...      % AKF-AR estimates for TPPSM
    ahl_kf_ar_tppsm, ...   % AHL-KF-AR estimates for TPPSM
    kf_std_tppsm, ...      % KF-std estimates for TPPSM
    akf_std_tppsm, ...     % AKF-std estimates for TPPSM
    ahl_kf_std_tppsm, ...  % AHL-KF-std estimates for TPPSM
    doppler_profile, ...
    cnr, ...
    seed, ...
    process_noise_variance, ...
    ar_order)
% plot_KFAR_estimates_comparison
%   Creates a uifigure containing a tab group with three tabs:
%       - "LOS Estimates"
%       - "Scintillation Estimates"
%       - "Joint Estimates"
%
%   In each tab, two subplots (one for the CSM scenario and one for the TPPSM
%   scenario) are arranged in a 2x1 tiled layout.
%
%   For AR-augmented estimates, the LOS component is in column 1 and the scintillation
%   component in column (length(doppler_profile)+1); their sum gives the joint estimate.
%   Standard KF (training_scint_model = 'none') only provides one column (used as the joint estimate).
%
%   A global title displaying simulation parameters (C/N0, Seed, Process Noise Variance,
%   AR Model Order) is added to the tab group. The entire uifigure is then saved as both
%   a PDF and a FIG file in a "figures" folder.
%
% Inputs:
%   time_vector        - [Nx1] vector with time stamps (sec).
%   los_phase          - [Nx1] vector with LOS phase (for detrending).
%   rx_sig_csm         - [Nx1] complex received signal for CSM.
%   kf_ar_csm          - [NxM] KF-AR estimates for CSM.
%   akf_ar_csm         - [NxM] AKF-AR estimates for CSM.
%   ahl_kf_ar_csm      - [NxM] AHL-KF-AR estimates for CSM.
%   kf_std_csm         - [NxM] KF-std estimates for CSM.
%   akf_std_csm        - [NxM] AKF-std estimates for CSM.
%   ahl_kf_std_csm     - [NxM] AHL-KF-std estimates for CSM.
%   rx_sig_tppsm       - [Nx1] complex received signal for TPPSM.
%   kf_ar_tppsm        - [NxM] KF-AR estimates for TPPSM.
%   akf_ar_tppsm       - [NxM] AKF-AR estimates for TPPSM.
%   ahl_kf_ar_tppsm    - [NxM] AHL-KF-AR estimates for TPPSM.
%   kf_std_tppsm       - [NxM] KF-std estimates for TPPSM.
%   akf_std_tppsm      - [NxM] AKF-std estimates for TPPSM.
%   ahl_kf_std_tppsm   - [NxM] AHL-KF-std estimates for TPPSM.
%   doppler_profile    - Vector used to determine the scintillation column (LOS is column 1).
%   cnr                - Carrier-to-noise ratio (dBHz).
%   seed               - Seed used in the simulation.
%   process_noise_variance - Process noise variance.
%   ar_order           - Order of the AR model.
%
% Author: [Your Name]
% Date: [Today's Date]

% Index for scintillation component (LOS is column 1)
scint_idx = length(doppler_profile) + 1;
% Global simulation parameters text
fprintf('Simulation Parameters: C/N0: %d dBHz; Seed: %d; Process Noise Variance: %g; AR Model Order: %d', ...
    cnr, seed, process_noise_variance, ar_order);

%% Create uifigure with tab group
fig = uifigure('Position',[50 50 1200 800]);
tabGroup = uitabgroup(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Create tabs for each type of estimate
tabLOS = uitab(tabGroup, 'Title', 'LOS Estimates');
tabScint = uitab(tabGroup, 'Title', 'Scintillation Estimates');
tabJoint = uitab(tabGroup, 'Title', 'Joint Estimates');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 1: LOS Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tloLOS = tiledlayout(tabLOS, 2, 1, 'TileSpacing', 'Compact','Padding','Compact');

% Subplot for CSM LOS estimates
ax1 = nexttile(tloLOS);
true_csm = unwrap(angle(rx_sig_csm)) - los_phase;
los_kf_ar = kf_ar_csm(:,1) - los_phase;
los_akf_ar = akf_ar_csm(:,1) - los_phase;
los_ahl_kf_ar = ahl_kf_ar_csm(:,1) - los_phase;
hold(ax1, 'on');
plot(ax1, time_vector, true_csm, 'k', 'LineWidth', 1.5);
plot(ax1, time_vector, los_kf_ar, 'r', 'LineWidth', 1.0);
plot(ax1, time_vector, los_akf_ar, 'b', 'LineWidth', 1.0);
plot(ax1, time_vector, los_ahl_kf_ar, 'g', 'LineWidth', 1.0);
hold(ax1, 'off');
grid(ax1, 'on');
xlabel(ax1, 'Time [sec]');
ylabel(ax1, 'Phase [rad]');
title(ax1, 'CSM - LOS Estimates');
legend(ax1, {'True Unwrapped Phase (LOS detrended)', 'KF-AR LOS','AKF-AR LOS','AHL-KF-AR LOS'}, ...
    'Location', 'best');

% Subplot for TPPSM LOS estimates
ax2 = nexttile(tloLOS);
true_tppsm = unwrap(angle(rx_sig_tppsm)) - los_phase;
los_kf_ar_tppsm = kf_ar_tppsm(:,1) - los_phase;
los_akf_ar_tppsm = akf_ar_tppsm(:,1) - los_phase;
los_ahl_kf_ar_tppsm = ahl_kf_ar_tppsm(:,1) - los_phase;
hold(ax2, 'on');
plot(ax2, time_vector, true_tppsm, 'k', 'LineWidth', 1.5);
plot(ax2, time_vector, los_kf_ar_tppsm, 'r', 'LineWidth', 1.0);
plot(ax2, time_vector, los_akf_ar_tppsm, 'b', 'LineWidth', 1.0);
plot(ax2, time_vector, los_ahl_kf_ar_tppsm, 'g', 'LineWidth', 1.0);
hold(ax2, 'off');
grid(ax2, 'on');
xlabel(ax2, 'Time [sec]');
ylabel(ax2, 'Phase [rad]');
title(ax2, 'TPPSM - LOS Estimates');
legend(ax2, {'True Unwrapped Phase (LOS detrended)','KF-AR LOS','AKF-AR LOS','AHL-KF-AR LOS'}, ...
    'Location', 'best');

save_tab(tabLOS, 'comparison_los');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 2: Scintillation Estimates (AR methods only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tloScint = tiledlayout(tabScint, 2, 1, 'TileSpacing', 'Compact','Padding','Compact');

% Subplot for CSM Scintillation estimates
ax3 = nexttile(tloScint);
true_wrapped_csm = angle(rx_sig_csm) - los_phase;
scint_kf_ar_csm = kf_ar_csm(:,scint_idx);
scint_akf_ar_csm = akf_ar_csm(:,scint_idx);
scint_ahl_kf_ar_csm = ahl_kf_ar_csm(:,scint_idx);
hold(ax3, 'on');
plot(ax3, time_vector, true_wrapped_csm, 'k', 'LineWidth', 1.5);
plot(ax3, time_vector, scint_kf_ar_csm, 'r', 'LineWidth', 1.0);
plot(ax3, time_vector, scint_akf_ar_csm, 'b', 'LineWidth', 1.0);
plot(ax3, time_vector, scint_ahl_kf_ar_csm, 'g', 'LineWidth', 1.0);
hold(ax3, 'off');
grid(ax3, 'on');
xlabel(ax3, 'Time [sec]');
ylabel(ax3, 'Phase [rad]');
title(ax3, 'CSM - Scintillation Estimates (AR only)');
legend(ax3, {'True Wrapped Phase (LOS detrended)', 'KF-AR scint','AKF-AR scint','AHL-KF-AR scint'}, 'Location', 'best');

% Subplot for TPPSM Scintillation estimates
ax4 = nexttile(tloScint);
true_wrapped_tppsm = angle(rx_sig_tppsm) - los_phase;
scint_kf_ar_tppsm = kf_ar_tppsm(:,scint_idx);
scint_akf_ar_tppsm = akf_ar_tppsm(:,scint_idx);
scint_ahl_kf_ar_tppsm = ahl_kf_ar_tppsm(:,scint_idx);
hold(ax4, 'on');
plot(ax4, time_vector, true_wrapped_tppsm, 'k', 'LineWidth', 1.5);
plot(ax4, time_vector, scint_kf_ar_tppsm, 'r', 'LineWidth', 1.0);
plot(ax4, time_vector, scint_akf_ar_tppsm, 'b', 'LineWidth', 1.0);
plot(ax4, time_vector, scint_ahl_kf_ar_tppsm, 'g', 'LineWidth', 1.0);
hold(ax4, 'off');
grid(ax4, 'on');
xlabel(ax4, 'Time [sec]');
ylabel(ax4, 'Phase [rad]');
title(ax4, 'TPPSM - Scintillation Estimates (AR only)');
legend(ax4, {'True Wrapped Phase (LOS detrended)', 'KF-AR scint','AKF-AR scint','AHL-KF-AR scint'}, 'Location', 'best');

save_tab(tabScint, 'comparison_scint');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 3: Joint Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tloJoint = tiledlayout(tabJoint, 2, 1, 'TileSpacing', 'Compact','Padding','Compact');

% Subplot for CSM Joint estimates
ax5 = nexttile(tloJoint);
joint_kf_ar_csm = (kf_ar_csm(:,1) - los_phase) + kf_ar_csm(:,scint_idx);
joint_akf_ar_csm = (akf_ar_csm(:,1) - los_phase) + akf_ar_csm(:,scint_idx);
joint_ahl_kf_ar_csm = (ahl_kf_ar_csm(:,1) - los_phase) + ahl_kf_ar_csm(:,scint_idx);
joint_kf_std_csm = kf_std_csm(:,1) - los_phase;
joint_akf_std_csm = akf_std_csm(:,1) - los_phase;
joint_ahl_kf_std_csm = ahl_kf_std_csm(:,1) - los_phase;
hold(ax5, 'on');
plot(ax5, time_vector, true_csm, 'k', 'LineWidth', 1.5);
plot(ax5, time_vector, joint_kf_ar_csm, 'r', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_akf_ar_csm, 'b', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_ahl_kf_ar_csm, 'g', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_kf_std_csm, '--r', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_akf_std_csm, '--b', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_ahl_kf_std_csm, '--g', 'LineWidth', 1.0);
hold(ax5, 'off');
grid(ax5, 'on');
xlabel(ax5, 'Time [sec]');
ylabel(ax5, 'Phase [rad]');
title(ax5, 'CSM - Joint Estimates');
legend(ax5, {'True Phase (LOS detrended)', 'KF-AR joint','AKF-AR joint','AHL-KF-AR joint',...
    'KF-std joint','AKF-std joint','AHL-KF-std joint'}, 'Location', 'best');

% Subplot for TPPSM Joint estimates
ax6 = nexttile(tloJoint);
joint_kf_ar_tppsm = (kf_ar_tppsm(:,1) - los_phase) + kf_ar_tppsm(:,scint_idx);
joint_akf_ar_tppsm = (akf_ar_tppsm(:,1) - los_phase) + akf_ar_tppsm(:,scint_idx);
joint_ahl_kf_ar_tppsm = (ahl_kf_ar_tppsm(:,1) - los_phase) + ahl_kf_ar_tppsm(:,scint_idx);
joint_kf_std_tppsm = kf_std_tppsm(:,1) - los_phase;
joint_akf_std_tppsm = akf_std_tppsm(:,1) - los_phase;
joint_ahl_kf_std_tppsm = ahl_kf_std_tppsm(:,1) - los_phase;
hold(ax6, 'on');
plot(ax6, time_vector, true_tppsm, 'k', 'LineWidth', 1.5);
plot(ax6, time_vector, joint_kf_ar_tppsm, 'r', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_akf_ar_tppsm, 'b', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_ahl_kf_ar_tppsm, 'g', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_kf_std_tppsm, '--r', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_akf_std_tppsm, '--b', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_ahl_kf_std_tppsm, '--g', 'LineWidth', 1.0);
hold(ax6, 'off');
grid(ax6, 'on');
xlabel(ax6, 'Time [sec]');
ylabel(ax6, 'Phase [rad]');
title(ax6, 'TPPSM - Joint Estimates');
legend(ax6, {'True Phase (LOS detrended)', 'KF-AR joint','AKF-AR joint','AHL-KF-AR joint',...
    'KF-std joint','AKF-std joint','AHL-KF-std joint'}, 'Location', 'best');

save_tab(tabJoint, 'comparison_joint');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helper function to save the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_tab(tabHandle, filename_prefix)
    % Get current folder and create figures folder if it doesn't exist
    current_folder = pwd;
    [baseFolder, ~, ~] = fileparts(current_folder);
    figures_dir = fullfile(baseFolder, 'figures');
    if ~exist(figures_dir, 'dir')
        mkdir(figures_dir);
    end

    % Force complete rendering
    drawnow;
    
    % Use exportgraphics to save the tab (container) as PDF (vector graphics)
    pdf_filename = fullfile(figures_dir, sprintf('%s.pdf', filename_prefix));
    exportgraphics(tabHandle, pdf_filename, 'ContentType', 'vector');
    
    % For saving as a FIG, we cannot directly save a tab.
    % Instead, we save the entire uifigure containing the tab.
    figHandle = ancestor(tabHandle, 'figure');
    fig_filename = fullfile(figures_dir, sprintf('%s.fig', filename_prefix));
    savefig(figHandle, fig_filename);
end

end

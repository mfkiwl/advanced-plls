function plot_time_series_full(...
    time_vector, ...
    los_phase, ...
    psi_csm, ...
    kf_ar_csm, ...         % KF-AR estimates (no adaptive update) for CSM
    akf_ar_csm, ...        % AKF-AR estimates (simplified, no HL) for CSM
    ahl_kf_ar_csm, ...     % AHL-KF-AR estimates (simplified, HL) for CSM
    kf_std_csm, ...        % KF-std estimates (no adaptive update) for CSM
    akf_std_csm, ...       % AKF-std estimates (simplified, no HL) for CSM
    ahl_kf_std_csm, ...    % AHL-KF-std estimates (simplified, HL) for CSM
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
    ar_order, ...
    refractive_phase_tppsm, ...
    diffractive_phase_tppsm, ...
    is_refractive_effects_removed_received_signal, ...
    is_refractive_effects_removed_training_data, ...
    is_save_figures)
% plot_time_series_full
%   Creates a uifigure containing a tab group with three tabs:
%       - "LOS Estimates"
%       - "Scintillation Estimates"
%       - "Joint Estimates"
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
%   psi_csm            - [Nx1] complex received signal for CSM.
%   kf_ar_csm          - [NxM] KF-AR estimates for CSM.
%   akf_ar_csm         - [NxM] AKF-AR estimates for CSM.
%   ahl_kf_ar_csm      - [NxM] AHL-KF-AR estimates for CSM.
%   kf_std_csm         - [NxM] KF-std estimates for CSM.
%   akf_std_csm        - [NxM] AKF-std estimates for CSM.
%   ahl_kf_std_csm     - [NxM] AHL-KF-std estimates for CSM.
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
%   refractive_phase_tppsm - [Nx1] refractive phase for TPPSM.
%   diffractive_phase_tppsm - [Nx1] diffractive phase for TPPSM.
%   is_refractive_effects_removed_received_signal - Boolean flag for received signal.
%   is_refractive_effects_removed_training_data   - Boolean flag for training data.
%   is_save_figures    - Flag to save the figures in pdf and .fig
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

if is_refractive_effects_removed_received_signal && ~is_refractive_effects_removed_training_data
    error("The situation when `is_refractive_effects_removed_received_signal` is true and `is_refractive_effects_removed_training_data` is not feasible.");
end

% Index for scintillation component (LOS is column 1)
scint_idx = length(doppler_profile) + 1;
% Global simulation parameters text
fprintf('Simulation Parameters: C/N0: %d dB-Hz; Seed: %d; Process Noise Variance: %g; AR Model Order: %d \n', ...
    cnr, seed, process_noise_variance, ar_order);

%% Create uifigure with tab group
fig = uifigure('Position',[50 50 1200 490]);
tabGroup = uitabgroup(fig, 'Units', 'normalized', 'Position', [0 0 1 1]);

% Create tabs for each type of estimate
tabLOS = uitab(tabGroup, 'Title', 'LOS Estimates');
tabScint = uitab(tabGroup, 'Title', 'Scintillation Estimates');
tabJoint = uitab(tabGroup, 'Title', 'Joint Estimates');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 1: LOS Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tloLOS = tiledlayout(tabLOS, 2, 1, 'TileSpacing', 'Compact','Padding','Compact');

%% Subplot for CSM LOS estimates
ax1 = nexttile(tloLOS);
% For CSM, plot both wrapped and unwrapped phases
true_wrapped_csm = angle(psi_csm);
true_unwrapped_csm = unwrap(angle(psi_csm));
los_detrended_kf_ar = kf_ar_csm(:,1) - los_phase;
los_detrended_akf_ar = akf_ar_csm(:,1) - los_phase;
los_detrended_ahl_kf_ar = ahl_kf_ar_csm(:,1) - los_phase;
hold(ax1, 'on');
plot(ax1, time_vector, true_wrapped_csm, 'k', 'LineWidth', 2.0);      % Wrapped CSM
plot(ax1, time_vector, true_unwrapped_csm, '--k', 'LineWidth', 2.0);   % Unwrapped CSM
plot(ax1, time_vector, los_detrended_kf_ar, 'r', 'LineWidth', 1.0);
plot(ax1, time_vector, los_detrended_akf_ar, 'b', 'LineWidth', 1.0);
plot(ax1, time_vector, los_detrended_ahl_kf_ar, 'g', 'LineWidth', 1.0);
hold(ax1, 'off');
grid(ax1, 'on');
xlabel(ax1, 'Time [sec]');
ylabel(ax1, 'Phase [rad]');
title(ax1, 'CSM - LOS Estimates (AR Only) - Without Zoom');
legend(ax1, {'Wrapped CSM','Unwrapped CSM','KF-AR LOS','AKF-AR LOS','AHL-KF-AR LOS'}, ...
    'Location', 'northwest');

%% Subplot for TPPSM LOS estimates
ax2 = nexttile(tloLOS);
los_detrended_kf_ar_tppsm = kf_ar_tppsm(:,1) - los_phase;
los_detrended_akf_ar_tppsm = akf_ar_tppsm(:,1) - los_phase;
los_detrended_ahl_kf_ar_tppsm = ahl_kf_ar_tppsm(:,1) - los_phase;
hold(ax2, 'on');
if is_refractive_effects_removed_received_signal && is_refractive_effects_removed_training_data
    % Case: Only diffractive phase is used (both flags true)
    tppsm_diffractive_wrapped = diffractive_phase_tppsm;
    tppsm_diffractive_unwrapped = unwrap(diffractive_phase_tppsm);
    plot(ax2, time_vector, tppsm_diffractive_wrapped, 'k', 'LineWidth', 2.0);
    plot(ax2, time_vector, tppsm_diffractive_unwrapped, '--k', 'LineWidth', 2.0);
    legend_entries = {'Diffractive Wrapped','Diffractive Unwrapped',...
        'KF-AR LOS','AKF-AR LOS','AHL-KF-AR LOS'};
    title_str = 'TPPSM - LOS Estimates (Diffractive Only) - Without Zoom';
else
    % Case: Received signal refractive effects NOT removed -> plot all components
    tppsm_refractive = refractive_phase_tppsm;
    tppsm_propagated_unwrapped = refractive_phase_tppsm + unwrap(diffractive_phase_tppsm);
    % Plot each reference signal with distinct linestyles/linewidths
    plot(ax2, time_vector, tppsm_refractive, ':k', 'LineWidth', 2.0);             % Refractive
    plot(ax2, time_vector, tppsm_propagated_unwrapped, '-.k', 'LineWidth', 3.0);      % Propagated unwrapped
    legend_entries = {'Refractive', 'Propagated Unwrapped', 'KF-AR LOS','AKF-AR LOS','AHL-KF-AR LOS'};
    title_str = 'TPPSM - LOS Estimates (Diffractive and Refractive) - Without Zoom';
end
% Now plot the Kalman filter LOS estimates on top
plot(ax2, time_vector, los_detrended_kf_ar_tppsm, 'r', 'LineWidth', 1.0);
plot(ax2, time_vector, los_detrended_akf_ar_tppsm, 'b', 'LineWidth', 1.0);
plot(ax2, time_vector, los_detrended_ahl_kf_ar_tppsm, 'g', 'LineWidth', 1.0);
hold(ax2, 'off');
grid(ax2, 'on');
xlabel(ax2, 'Time [sec]');
ylabel(ax2, 'Phase [rad]');
title(ax2, title_str);
legend(ax2, legend_entries, 'Location', 'best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 2: Scintillation Estimates (AR methods only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tloScint = tiledlayout(tabScint, 2, 1, 'TileSpacing', 'Compact','Padding','Compact');

%% Subplot for CSM Scintillation estimates
ax3 = nexttile(tloScint);
% For CSM, plot both wrapped and unwrapped phases as reference
true_wrapped_csm = angle(psi_csm);
true_unwrapped_csm = unwrap(true_wrapped_csm);
scint_kf_ar_csm = kf_ar_csm(:,scint_idx);
scint_akf_ar_csm = akf_ar_csm(:,scint_idx);
scint_ahl_kf_ar_csm = ahl_kf_ar_csm(:,scint_idx);
hold(ax3, 'on');
plot(ax3, time_vector, true_wrapped_csm, 'k', 'LineWidth', 2.0);
plot(ax3, time_vector, true_unwrapped_csm, '--k', 'LineWidth', 2.0);
plot(ax3, time_vector, scint_kf_ar_csm, 'r', 'LineWidth', 1.0);
plot(ax3, time_vector, scint_akf_ar_csm, 'b', 'LineWidth', 1.0);
plot(ax3, time_vector, scint_ahl_kf_ar_csm, 'g', 'LineWidth', 1.0);
hold(ax3, 'off');
grid(ax3, 'on');
xlabel(ax3, 'Time [sec]');
ylabel(ax3, 'Phase [rad]');
title(ax3, 'CSM - Scintillation Estimates (AR Only) - Without Zoom');
legend(ax3, {'Wrapped CSM','Unwrapped CSM','KF-AR Scint','AKF-AR Scint','AHL-KF-AR Scint'}, 'Location', 'southwest');

%% Subplot for TPPSM Scintillation estimates
ax4 = nexttile(tloScint);
hold(ax4, 'on');
if is_refractive_effects_removed_received_signal && is_refractive_effects_removed_training_data
    % Case: Only diffractive phase is used for scintillation reference
    tppsm_diffractive_wrapped = diffractive_phase_tppsm;
    tppsm_diffractive_unwrapped = unwrap(diffractive_phase_tppsm);
    plot(ax4, time_vector, tppsm_diffractive_wrapped, 'k', 'LineWidth', 2.0);
    plot(ax4, time_vector, tppsm_diffractive_unwrapped, '--k', 'LineWidth', 2.0);
    legend_entries_joint = {'Diffractive Wrapped','Diffractive Unwrapped',...
        'KF-AR Scint','AKF-AR Scint','AHL-KF-AR Scint'};
    title_str_joint = 'TPPSM - Scintillation Estimates (Diffractive Only) - Without Zoom';
else
    tppsm_propagated_wrapped = wrapToPi(refractive_phase_tppsm + diffractive_phase_tppsm);
    tppsm_propagated_unwrapped = refractive_phase_tppsm + unwrap(diffractive_phase_tppsm);
    tppsm_diffractive_unwrapped = unwrap(diffractive_phase_tppsm);
    plot(ax4, time_vector, diffractive_phase_tppsm, '-k', 'LineWidth', 2.0);
    plot(ax4, time_vector, tppsm_diffractive_unwrapped, ':k', 'LineWidth', 2.0);
    plot(ax4, time_vector, tppsm_propagated_wrapped, '--k', 'LineWidth', 2.0);
    plot(ax4, time_vector, tppsm_propagated_unwrapped, '-.k', 'LineWidth', 2.0);
    legend_entries_joint = {'Diffractive Wrapped', 'Diffractive Unwrapped','Propagated Wrapped','Propagated Unwrapped','KF-AR Scint','AKF-AR Scint','AHL-KF-AR Scint'};
    title_str_joint = 'TPPSM - Scintillation Estimates (Diffractive and Refractive) - Without Zoom';
end
% Now overlay the Kalman filter scintillation estimates
scint_kf_ar_tppsm = kf_ar_tppsm(:,scint_idx);
scint_akf_ar_tppsm = akf_ar_tppsm(:,scint_idx);
scint_ahl_kf_ar_tppsm = ahl_kf_ar_tppsm(:,scint_idx);
plot(ax4, time_vector, scint_kf_ar_tppsm, 'r', 'LineWidth', 1.0);
plot(ax4, time_vector, scint_akf_ar_tppsm, 'b', 'LineWidth', 1.0);
plot(ax4, time_vector, scint_ahl_kf_ar_tppsm, 'g', 'LineWidth', 1.0);
hold(ax4, 'off');
grid(ax4, 'on');
xlabel(ax4, 'Time [sec]');
ylabel(ax4, 'Phase [rad]');
title(ax4, title_str_joint);
legend(ax4, legend_entries_joint, 'Location', 'southwest');

if is_save_figures
    save_tab(tabScint, 'comparison_scint');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tab 3: Joint Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tloJoint = tiledlayout(tabJoint, 2, 1, 'TileSpacing', 'Compact','Padding','Compact');

%% Subplot for CSM Joint estimates
ax5 = nexttile(tloJoint);
joint_detrended_kf_ar_csm = (kf_ar_csm(:,1) - los_phase) + kf_ar_csm(:,scint_idx);
joint_detrended_akf_ar_csm = (akf_ar_csm(:,1) - los_phase) + akf_ar_csm(:,scint_idx);
joint_detrended_ahl_kf_ar_csm = (ahl_kf_ar_csm(:,1) - los_phase) + ahl_kf_ar_csm(:,scint_idx);
joint_detrended_kf_std_csm = kf_std_csm(:,1) - los_phase;
joint_detrended_akf_std_csm = akf_std_csm(:,1) - los_phase;
joint_detrended_ahl_kf_std_csm = ahl_kf_std_csm(:,1) - los_phase;
hold(ax5, 'on');
plot(ax5, time_vector, true_wrapped_csm, '--k', 'LineWidth', 2.0);
plot(ax5, time_vector, true_unwrapped_csm, 'k', 'LineWidth', 2.0);
plot(ax5, time_vector, joint_detrended_kf_ar_csm, 'r', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_detrended_akf_ar_csm, 'b', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_detrended_ahl_kf_ar_csm, 'g', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_detrended_kf_std_csm, '--r', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_detrended_akf_std_csm, '--b', 'LineWidth', 1.0);
plot(ax5, time_vector, joint_detrended_ahl_kf_std_csm, '--g', 'LineWidth', 1.0);
hold(ax5, 'off');
grid(ax5, 'on');
xlabel(ax5, 'Time [sec]');
ylabel(ax5, 'Phase [rad]');
title(ax5, 'CSM - Joint Estimates - Without Zoom');
legend(ax5, {'True Wrapped Phase','True Unwrapped Phase','KF-AR joint','AKF-AR joint','AHL-KF-AR joint',...
    'KF-std joint','AKF-std joint','AHL-KF-std joint'}, 'Location', 'northwest');

%% Subplot for TPPSM Joint estimates
ax6 = nexttile(tloJoint);
joint_detrended_kf_ar_tppsm = (kf_ar_tppsm(:,1) - los_phase) + kf_ar_tppsm(:,scint_idx);
joint_detrended_akf_ar_tppsm = (akf_ar_tppsm(:,1) - los_phase) + kf_ar_tppsm(:,scint_idx);
joint_detrended_ahl_kf_ar_tppsm = (ahl_kf_ar_tppsm(:,1) - los_phase) + ahl_kf_ar_tppsm(:,scint_idx);
joint_detrended_kf_std_tppsm = kf_std_tppsm(:,1) - los_phase;
joint_detrended_akf_std_tppsm = akf_std_tppsm(:,1) - los_phase;
joint_detrended_ahl_kf_std_tppsm = ahl_kf_std_tppsm(:,1) - los_phase;
if is_refractive_effects_removed_received_signal && is_refractive_effects_removed_training_data
    % Case: Only diffractive phase is used for scintillation reference
    tppsm_diffractive_wrapped = diffractive_phase_tppsm;
    tppsm_diffractive_unwrapped = unwrap(diffractive_phase_tppsm);
    plot(ax6, time_vector, tppsm_diffractive_wrapped, 'k', 'LineWidth', 2.0);
    plot(ax6, time_vector, tppsm_diffractive_unwrapped, '--k', 'LineWidth', 2.0);
    legend_entries_joint = {'Diffractive Wrapped','Diffractive Unwrapped',...
        'KF-AR Scint','AKF-AR Scint','AHL-KF-AR Scint'};
    title_str_joint = 'TPPSM - Scintillation Estimates (Diffractive Only) - Without Zoom';
else
    % Case: Refractive phase is present on the received signal
    tppsm_propagated_unwrapped = refractive_phase_tppsm + unwrap(diffractive_phase_tppsm);
    plot(ax6, time_vector, tppsm_propagated_unwrapped, ':k', 'LineWidth', 3.0);
    legend_entries_joint = {'Propagated Unwrapped','KF-AR Scint','AKF-AR Scint','AHL-KF-AR Scint'};
    title_str_joint = 'TPPSM - Scintillation Estimates (Diffractive and Refractive) - Without Zoom';
end
hold(ax6, 'on');
plot(ax6, time_vector, joint_detrended_kf_ar_tppsm, 'r', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_detrended_akf_ar_tppsm, 'b', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_detrended_ahl_kf_ar_tppsm, 'g', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_detrended_kf_std_tppsm, '--r', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_detrended_akf_std_tppsm, '--b', 'LineWidth', 1.0);
plot(ax6, time_vector, joint_detrended_ahl_kf_std_tppsm, '--g', 'LineWidth', 1.0);
hold(ax6, 'off');
grid(ax6, 'on');
xlabel(ax6, 'Time [sec]');
ylabel(ax6, 'Phase [rad]');
title(ax6, title_str_joint);
legend(ax6, legend_entries_joint, 'Location', 'best');
if is_save_figures
    save_tab(tabJoint, 'comparison_joint');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local helper function to save the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_tab(tabHandle, filename_prefix)
    % Get current folder and create figures folder if it doesn't exist
    current_folder = pwd;
    [baseFolder, ~, ~] = fileparts(current_folder);
    figures_dir = fullfile(baseFolder, 'figures', 'time_series_comparison');
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

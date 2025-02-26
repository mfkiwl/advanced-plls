function plot_rmse_heatmaps(rmse_csm, rmse_tppsm, process_noise_variance_array, phase_type, algorithm_label, directory_to_save)
% plot_rmse_heatmaps
% Plots side-by-side 2x1 heatmaps comparing the RMSE distribution versus process
% noise variance for two scenarios (CSM and TPPSM) for a given phase type and algorithm.
%
% This function builds a 2D histogram (using fixed RMSE bins) for the RMSE values at each
% noise level and displays the resulting heatmap. Overlaid on the heatmap are white lines:
%   - A solid line representing the mean RMSE at each noise level.
%   - Dashed lines representing the 5th and 95th percentile (confidence bounds).
%
% Inputs:
%   rmse_csm       - [N_mc x N_noise] RMSE values for the CSM scenario.
%   rmse_tppsm     - [N_mc x N_noise] RMSE values for the TPPSM scenario.
%   noise_vals     - Vector of process noise variance values (length = N_noise).
%   phase_type     - String describing the phase type (e.g., 'los', 'scint', 'joint').
%   algorithm_label- String describing the algorithm variant (e.g., 'kf_ar', 'akf_ar', 'ahl_kf_ar', etc.).
%
% Example:
%   plot_rmse_heatmaps(results.csm.los.kf_ar, results.tppsm.los.kf_ar, pnv_array, 'los', 'kf_ar');
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

% Determine RMSE bins based on combined data (for consistency)
all_data = [rmse_csm(:); rmse_tppsm(:)];
min_rmse = min(all_data);
max_rmse = max(all_data);
num_bins = 50;
edges_rmse = linspace(min_rmse, max_rmse, num_bins+1);
centers_rmse = (edges_rmse(1:end-1) + edges_rmse(2:end)) / 2;

% Create figure and tiled layout
fig = figure('Position',[50,50,500,600]);
tiled_layout = tiledlayout(2, 1, 'TileSpacing', 'Compact', 'Padding', 'Compact');

%% Panel 1: CSM
nexttile(tiled_layout);
heatmap_data = zeros(num_bins, length(process_noise_variance_array));
for j = 1:length(process_noise_variance_array)
    counts = histcounts(rmse_csm(:,j), edges_rmse);
    heatmap_data(:,j) = counts(:);
end
imagesc(process_noise_variance_array, centers_rmse, heatmap_data);
set(gca, 'XScale', 'log'); % Use log-scale for noise axis
colorbar;
xlabel('Process Noise Variance');
ylabel('RMSE');
title(sprintf('%s %s: CSM', algorithm_label, phase_type));
hold on;
% Compute summary statistics per noise level
mean_rmse = mean(rmse_csm, 1);
lower_bound = prctile(rmse_csm, 5, 1);
upper_bound = prctile(rmse_csm, 95, 1);
plot(process_noise_variance_array, mean_rmse, 'w-', 'LineWidth', 2);   % white solid line for mean
plot(process_noise_variance_array, lower_bound, 'w--', 'LineWidth', 2);  % white dashed line for lower bound
plot(process_noise_variance_array, upper_bound, 'w--', 'LineWidth', 2);  % white dashed line for upper bound
hold off;

%% Panel 2: TPPSM
nexttile(tiled_layout);
heatmap_data = zeros(num_bins, length(process_noise_variance_array));
for j = 1:length(process_noise_variance_array)
    counts = histcounts(rmse_tppsm(:,j), edges_rmse);
    heatmap_data(:,j) = counts(:);
end
imagesc(process_noise_variance_array, centers_rmse, heatmap_data);
set(gca, 'XScale', 'log');
colorbar;
xlabel('Process Noise Variance');
ylabel('RMSE');
title(sprintf('%s %s: TPPSM', algorithm_label, phase_type));
hold on;
mean_rmse = mean(rmse_tppsm, 1);
lower_bound = prctile(rmse_tppsm, 5, 1);
upper_bound = prctile(rmse_tppsm, 95, 1);
plot(process_noise_variance_array, mean_rmse, 'w-', 'LineWidth', 2);
plot(process_noise_variance_array, lower_bound, 'w--', 'LineWidth', 2);
plot(process_noise_variance_array, upper_bound, 'w--', 'LineWidth', 2);
hold off;

sgtitle(sprintf('%s %s RMSE vs. Process Noise Variance', algorithm_label, phase_type));

% Adjust figure paper size
set(fig, 'PaperUnits', 'centimeters');
fig_pos = get(fig, 'Position');
fig_width_cm = fig_pos(3) * 2.54 / 96;
fig_height_cm = fig_pos(4) * 2.54 / 96;
set(fig, 'PaperSize', [fig_width_cm fig_height_cm]);  
set(fig, 'PaperPosition', [0 0 fig_width_cm fig_height_cm]);

% Save as PDF
pdf_filename = sprintf('%s_%s_RMSE_hist_pnv.pdf', algorithm_label, phase_type);
print(fig, fullfile(directory_to_save, pdf_filename), '-dpdf', '-vector');

disp(['Saved: ', fullfile(directory_to_save, pdf_filename)]);
close(fig);
end

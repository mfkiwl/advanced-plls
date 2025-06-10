% Script: plot_all_results.m
%
% Description:
%   This script loads the RMSE results obtained from the Monte Carlo
%   simulation with process noise sweep (stored in a 'results' folder in the
%   parent directory) and generates side-by-side 2x1 heatmaps comparing the RMSE
%   distributions versus process noise variance for different phase types and
%   algorithm variants.
%
%   The following comparisons are made:
%
%   1. LOS RMSE:
%      1.1 - Compare KF-AR los phase RMSE in CSM and TPPSM.
%      1.2 - Compare AKF-AR los phase RMSE in CSM and TPPSM.
%      1.3 - Compare AHL-KF-AR los phase RMSE in CSM and TPPSM.
%
%   2. Scintillation RMSE (AR methods only):
%      2.1 - Compare KF-AR scint phase RMSE in CSM and TPPSM.
%      2.2 - Compare AKF-AR scint phase RMSE in CSM and TPPSM.
%      2.3 - Compare AHL-KF-AR scint phase RMSE in CSM and TPPSM.
%
%   3. Joint RMSE:
%      3.1 - Compare KF-AR joint phase RMSE in CSM and TPPSM.
%      3.2 - Compare AKF-AR joint phase RMSE in CSM and TPPSM.
%      3.3 - Compare AHL-KF-AR joint phase RMSE in CSM and TPPSM.
%      3.4 - Compare KF-std joint phase RMSE in CSM and TPPSM.
%      3.5 - Compare AKF-std joint phase RMSE in CSM and TPPSM.
%      3.6 - Compare AHL-KF-std joint phase RMSE in CSM and TPPSM.
%
%   Each plot is produced as a 2x1 tiled layout (top panel for CSM and bottom panel for TPPSM).
%   Over each heatmap, white lines indicate the mean RMSE and the 5th and 95th percentile
%   (confidence margins).
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clear; clc;

%% Load Results Structure
script_folder = fileparts(mfilename('fullpath'));
results_folder = fullfile(script_folder, 'results');
results_file = fullfile(results_folder, 'monte_carlo_with_process_noise_sweep_results.mat');
libs_folder = genpath(fullfile('..','..','libs'));

addpath(libs_folder);

if exist(results_file, 'file')
    loaded_data = load(results_file, 'results');
    results = loaded_data.results;
else
    error('plot_all_results:resultsNotFound', 'Results file not found: %s', results_file);
end

directory_to_save = fullfile(fileparts(fileparts(mfilename('fullpath'))),'figures','RMSE_hist_vs_pnv');

if ~exist(directory_to_save, 'dir')
    mkdir(directory_to_save);
end

% Determine the number of noise levels from the results structure
process_noise_variance_array = results.pnv_array;

%% Plotting

% --- 1. LOS RMSE Heatmaps ---
% 1.1: KF-AR los phase RMSE
plot_rmse_heatmaps(results.csm.los.kf_ar, results.tppsm.los.kf_ar, process_noise_variance_array, 'line-of-sight phase', 'KF-AR', directory_to_save);
% 1.2: AKF-AR los phase RMSE
plot_rmse_heatmaps(results.csm.los.akf_ar, results.tppsm.los.akf_ar, process_noise_variance_array, 'line-of-sight phase', 'AKF-AR', directory_to_save);
% 1.3: AHL-KF-AR los phase RMSE
plot_rmse_heatmaps(results.csm.los.ahl_kf_ar, results.tppsm.los.ahl_kf_ar, process_noise_variance_array, 'line-of-sight phase', 'AHL-KF-AR', directory_to_save);

% --- 2. Scintillation RMSE Heatmaps (AR estimates only) ---
% 2.1: KF-AR scint phase RMSE
plot_rmse_heatmaps(results.csm.scint.kf_ar, results.tppsm.scint.kf_ar, process_noise_variance_array, 'scintillation phase', 'KF-AR', directory_to_save);
% 2.2: AKF-AR scint phase RMSE
plot_rmse_heatmaps(results.csm.scint.akf_ar, results.tppsm.scint.akf_ar, process_noise_variance_array, 'scintillation phase', 'AKF-AR', directory_to_save);
% 2.3: AHL-KF-AR scint phase RMSE
plot_rmse_heatmaps(results.csm.scint.ahl_kf_ar, results.tppsm.scint.ahl_kf_ar, process_noise_variance_array, 'scintillation phase', 'AHL-KF-AR', directory_to_save);

% --- 3. Joint RMSE Heatmaps ---
% 3.1: KF-AR joint phase RMSE
plot_rmse_heatmaps(results.csm.joint.kf_ar, results.tppsm.joint.kf_ar, process_noise_variance_array, 'joint (LOS + scint) phase', 'KF-AR', directory_to_save);
% 3.2: AKF-AR joint phase RMSE
plot_rmse_heatmaps(results.csm.joint.akf_ar, results.tppsm.joint.akf_ar, process_noise_variance_array, 'joint (LOS + scint) phase', 'AKF-AR', directory_to_save);
% 3.3: AHL-KF-AR joint phase RMSE
plot_rmse_heatmaps(results.csm.joint.ahl_kf_ar, results.tppsm.joint.ahl_kf_ar, process_noise_variance_array, 'joint (LOS + scint) phase', 'AHL-KF-AR', directory_to_save);
% 3.4: KF-std joint phase RMSE
plot_rmse_heatmaps(results.csm.joint.kf_std, results.tppsm.joint.kf_std, process_noise_variance_array, 'joint (LOS + scint) phase', 'KF-Std', directory_to_save);
% 3.5: AKF-std joint phase RMSE
plot_rmse_heatmaps(results.csm.joint.akf_std, results.tppsm.joint.akf_std, process_noise_variance_array, 'joint (LOS + scint) phase', 'AKF-Std', directory_to_save);
% 3.6: AHL-KF-std joint phase RMSE
plot_rmse_heatmaps(results.csm.joint.ahl_kf_std, results.tppsm.joint.ahl_kf_std, process_noise_variance_array, 'joint (LOS + scint) phase', 'AHL-KF-Std', directory_to_save);

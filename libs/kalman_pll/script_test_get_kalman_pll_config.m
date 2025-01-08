% Clear workspace and command window
clear all; clc;

addpath(genpath(fullfile(pwd,'..')));

% Define discrete Wiener model configuration
discrete_wiener_model_config = {1,3,0.01,[0,0,1],1}; % Example: L, M, sampling_interval, sigma_array, delta_array
% Define scintillation training data configuration
scintillation_training_data_config = {0.8, 0.7, 300, 0.01}; % Example: S4, tau0, simulation_time, sampling_interval

% Set VAR model parameters
var_minimum_order = 1; % Minimum VAR model order
var_maximum_order = 6; % Maximum VAR model order

% Define C/N0 array
C_over_N0_array_dBHz = 35; % Example values in dB-Hz

% Choose scintillation model: 'CSM', 'MFPSM', or 'none'
training_scint_model = 'none';

% Set flag to remove refractive effects for MFPSM (not applicable for CSM or none)
is_refractive_effects_removed = true;

% Use cached models (set to true to load from cache, false to compute fresh)
is_use_cached_settings = false;

% Combine all configurations into a single struct
config = struct( ...
    'discrete_wiener_model_config', {discrete_wiener_model_config}, ...
    'scintillation_training_data_config', {scintillation_training_data_config}, ...
    'var_minimum_order', var_minimum_order, ...
    'var_maximum_order', var_maximum_order, ...
    'C_over_N0_array_dBHz', C_over_N0_array_dBHz, ...
    'training_scint_model', training_scint_model, ...
    'is_refractive_effects_removed', is_refractive_effects_removed, ...
    'is_use_cached_settings', is_use_cached_settings ...
);

kalman_pll_config = get_kalman_pll_config(config);


% Clear workspace and command window
clear all; clc;

addpath(genpath(fullfile(pwd,'..')));

% Set VAR model parameters
var_minimum_order = 1; % Minimum VAR model order
var_maximum_order = 6; % Maximum VAR model order

% Define C/N0 array
C_over_N0_array_dBHz = 45; % Example values in dB-Hz

% Choose scintillation model: 'CSM', 'MFPSM', or 'none'
training_scint_model = 'none';

% Set the initial states uniform distributions boundaries in a cell array
initial_states_distributions_boundaries = {[-pi,pi],[-5,5],[-0.0001,0.0001]};

% Set the Doppler profile used to simulate the real line-of-sight phase
% dynamics on the function `get_los_phase`
% TODO: Make a validation that assures that
% `initial_states_distributions_boundaries` and `real_doppler_profile` have
% the same amount of elements.
real_doppler_profile = [0, 1000, 0.94];

% Set flag to remove refractive effects for MFPSM (not applicable for CSM or none)
is_refractive_effects_removed = false;

% Set flag to use cached settings (set to true to load from cache, false to compute fresh)
is_use_cached_settings = false;

% Set flag to generate initial estimates
is_generate_random_initial_estimates = false;

S4 = 0.8;
tau0 = 0.7;
simulation_time = 300;
var_model_training_simulation_time = 900;
settling_time = 50;
scint_model = 'none';
sampling_interval = 0.01;

L = 1;
M = 3;
sigma_array = [0,0,1e-2];
delta_array = 1;

% Define discrete Wiener model configuration
discrete_wiener_model_config = {L,M,sampling_interval,sigma_array,delta_array}; % Example: L, M, sampling_interval, sigma_array, delta_array
% Define scintillation training data configuration
scintillation_training_data_config = {S4, tau0, simulation_time, sampling_interval}; % Example: S4, tau0, simulation_time, sampling_interval

% Combine all configurations into a single struct
config = struct( ...
    'discrete_wiener_model_config', {discrete_wiener_model_config}, ...
    'scintillation_training_data_config', {scintillation_training_data_config}, ...
    'var_minimum_order', var_minimum_order, ...
    'var_maximum_order', var_maximum_order, ...
    'C_over_N0_array_dBHz', C_over_N0_array_dBHz, ...
    'training_scint_model', training_scint_model, ...
    'initial_states_distributions_boundaries', {initial_states_distributions_boundaries}, ...
    'real_doppler_profile', real_doppler_profile, ...
    'is_refractive_effects_removed', is_refractive_effects_removed, ...
    'is_use_cached_settings', is_use_cached_settings, ...
    'is_generate_random_initial_estimates', is_generate_random_initial_estimates ...
);

time = 0.01:0.01:simulation_time;

[received_signal, los_phase, psi, ps_realization] = ...
    get_received_signal(C_over_N0_array_dBHz(1),S4,tau0,simulation_time, ...
    settling_time,scint_model, real_doppler_profile, is_refractive_effects_removed);

[kalman_pll_config,initial_estimates] = get_kalman_pll_config(config);

[state_estimates, error_covariance_estimates] = ...
    get_kalman_pll_estimates(received_signal,kalman_pll_config,initial_estimates,training_scint_model);

plot(time,los_phase - state_estimates(:,1));
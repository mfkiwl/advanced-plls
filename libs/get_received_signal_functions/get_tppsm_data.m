function [psi_tppsm, ps_realization, general_params, irr_params_set, seed] = get_tppsm_data(scenario, varargin)
% get_tppsm_data
% Generates Two-Component Power Law Model (TPPSM) realizations using preset 
% irregularity parameters for a specified scintillation scenario.
%
% TPPSM stands for Two-Component Power Law Model, and this implementation is 
% focused on single-frequency (e.g., L1) carrier phase tracking. This function 
% does not perform any extrapolation of the irregularity parameters; it uses 
% the preset parameters directly. The computed rhof/veff ratio for the L1 frequency
% is passed directly to get_scintillation_time_series.
%
% Syntax:
%   [psi_tppsm, ps_realization, general_params, irr_params_set, seed] = ...
%       get_tppsm_data(scenario, 'simulation_time', simulation_time, ...
%                      'sampling_interval', sampling_interval, ...
%                      'general_params', general_params, ...
%                      'irr_params_set', irr_params_set, ...
%                      'seed', seed)
%
% Inputs:
%   scenario         - A string specifying the scintillation scenario 
%                      ('Weak', 'Moderate', or 'Severe')
%
% Optional Name-Value Pair Inputs:
%   'simulation_time'  - Duration of the simulation in seconds (default: 300 s)
%   'sampling_interval'- Sampling time in seconds (default: 0.01 s)
%   'general_params' - Structure containing simulation parameters. If empty,
%                      the function will call get_general_parameters() internally.
%   'irr_params_set' - Structure containing preset irregularity parameters.
%                      If empty, the function will call get_irregularity_parameters().
%   'seed'           - Seed for random number generation (default: 1)
%
% Outputs:
%   psi_tppsm        - Scintillation field time series (complex)
%   ps_realization   - Phase screen (detrended phase) realization
%   general_params   - Simulation parameters used for this run.
%   irr_params_set   - Preset irregularity parameters used.
%   seed             - Seed value used.
%
% Example:
%   [psi_tppsm, ps_realization] = get_tppsm_data('Moderate', 'seed', 42);
%
% Notes:
%   - This implementation is for single-frequency carrier phase tracking.
%   - The function does not extrapolate irregularity parameters; instead, it uses
%     the preset parameters directly.
%   - The computed rhof/veff ratio for the L1 frequency is passed directly to
%     get_scintillation_time_series.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

% Input validation using inputParser
p = inputParser;
addRequired(p, 'scenario', @(x) ischar(x) || isstring(x));
addParameter(p, 'simulation_time', 300, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p, 'sampling_interval', 0.01, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(p, 'general_params', [], @(x) isstruct(x));
addParameter(p, 'irr_params_set', [], @(x) isstruct(x));
addParameter(p, 'seed', 1, @(x) isnumeric(x) && isscalar(x));
parse(p, scenario, varargin{:});

scenario = validatestring(p.Results.scenario, {'Weak', 'Moderate', 'Severe'}, mfilename, 'scenario');
simulation_time = p.Results.simulation_time;
sampling_interval = p.Results.sampling_interval;
general_params = p.Results.general_params;
irr_params_set = p.Results.irr_params_set;
seed = p.Results.seed;

if simulation_time < sampling_interval
    error('get_tppsm_data:simulationTimeSmallerThanSamplingInterval', ...
          'The simulation_time (%g s) is smaller than the sampling_interval (%g s).', ...
          simulation_time, sampling_interval);
end

% Compute number of samples
num_samples_exact = simulation_time / sampling_interval;
num_samples_rounded = round(num_samples_exact);
if abs(num_samples_exact - num_samples_rounded) > eps
    warning('get_tppsm_data:NonIntegerRatioSamples', ...
            'simulation_time / sampling_interval is not an integer. Rounded from %.5g to %d samples.', ...
            num_samples_exact, num_samples_rounded);
end

% Get or set simulation parameters
if isempty(general_params)
    general_params = get_general_parameters();
end
general_params.simulation_time = simulation_time;
general_params.dt = sampling_interval;

if isempty(irr_params_set)
    irr_params_set = get_irregularity_parameters();
end

% Ensure the selected scenario exists in the irregularity parameters
if ~isfield(irr_params_set, scenario)
    error('get_tppsm_data:InvalidScenario', ...
          'The scenario "%s" is not available in the preset irregularity parameters.', scenario);
end

% Use the preset irregularity parameters directly (no extrapolation)
irr_params = irr_params_set.(scenario);

% Compute the rhof/veff ratio for the L1 frequency (assuming L1 is primary)
rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_params);

% Set the random seed for reproducibility
rng(seed);

% Generate scintillation time series using the refactored interface.
% get_scintillation_time_series returns:
% [scint_field, norm_phase_sdf, detrended_phase, mu, doppler_frequency]
[scint_field, ~, detrended_phase, ~, ~] = get_scintillation_time_series(...
    general_params, irr_params, rhof_veff_ratio_L1, seed);

% Trim outputs to match the expected number of samples
psi_tppsm = scint_field(1:num_samples_rounded).';
ps_realization = detrended_phase(1:num_samples_rounded).';
end

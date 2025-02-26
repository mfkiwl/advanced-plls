function [psi_tppsm, ps_realization, general_params, irr_params_set, seed] = get_tppsm_data(scenario, varargin)
% get_tppsm_data
% Generates Two-Component Power Law Model (TPPSM) realizations using preset 
% irregularity parameters for a specified scintillation scenario.
%
% TPPSM stands for Two-Component Power Law Model. This implementation is 
% focused on single-frequency (e.g., L1) carrier phase tracking. The function 
% uses preset irregularity parameters (or those provided via the name-value pairs)
% and passes the computed or provided rhof/veff ratio for the L1 frequency directly
% to get_scintillation_time_series.
%
% Syntax:
%   [psi_tppsm, ps_realization, general_params, irr_params_set, seed] = ...
%       get_tppsm_data(scenario, 'simulation_time', simulation_time, ...
%                      'sampling_interval', sampling_interval, ...
%                      'general_params', general_params, ...
%                      'irr_params_set', irr_params_set, ...
%                      'seed', seed, ...
%                      'rhof_veff_ratio_L1', rhof_veff_ratio_L1)
%
% Inputs:
%   scenario         - A string specifying the scintillation scenario 
%                      ('Weak', 'Moderate', or 'Severe').
%
% Optional Name-Value Pair Inputs:
%   'simulation_time'     - Duration of the simulation in seconds (default: 300 s).
%   'sampling_interval'   - Sampling time in seconds (default: 0.01 s).
%   'general_params'      - Structure containing simulation parameters. If empty,
%                           the function will call get_general_parameters() internally.
%   'irr_params_set'      - Structure containing preset irregularity parameters.
%                           If empty, the function will call get_irregularity_parameters().
%   'seed'                - Seed for random number generation (default: 1).
%   'rhof_veff_ratio_L1'  - (Optional) A numeric value for the computed 
%                           rhof/veff ratio for L1. If provided (and nonempty), 
%                           the function will skip its internal computation.
%
% Outputs:
%   psi_tppsm        - Scintillation field time series (complex), as a row vector.
%   ps_realization   - Detrended phase screen realization (row vector).
%   general_params   - Simulation parameters used for this run.
%   irr_params_set   - Preset irregularity parameters used.
%   seed             - Seed value used.
%
% Example:
%   [psi_tppsm, ps_realization] = get_tppsm_data('Moderate', 'seed', 42, 'rhof_veff_ratio_L1', 0.35);
%
% Notes:
%   - This implementation is for single-frequency carrier phase tracking.
%   - The function does not perform any extrapolation of irregularity parameters;
%     instead, it uses preset parameters directly.
%   - If 'rhof_veff_ratio_L1' is provided, it is used directly; otherwise, the
%     function computes it by calling get_rhof_veff_ratio(general_params).
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

% Input validation using inputParser
p = inputParser;
addRequired(p, 'scenario', @(x) ischar(x) || isstring(x));
addParameter(p, 'simulation_time', 300, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'positive'}));
addParameter(p, 'sampling_interval', 0.01, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real', 'positive'}));
addParameter(p, 'general_params', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'irr_params_set', [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'seed', 1, @(x) validateattributes(x, {'numeric'}, {'scalar', 'real'}));
addParameter(p, 'rhof_veff_ratio_L1', [], @(x) isempty(x) || isnumeric(x));
parse(p, scenario, varargin{:});

scenario = validatestring(p.Results.scenario, {'Weak', 'Moderate', 'Severe'}, mfilename, 'scenario');
simulation_time = p.Results.simulation_time;
sampling_interval = p.Results.sampling_interval;
general_params = p.Results.general_params;
irr_params_set = p.Results.irr_params_set;
seed = p.Results.seed;
rhof_veff_ratio_L1 = p.Results.rhof_veff_ratio_L1;

if simulation_time < sampling_interval
    error('get_tppsm_data:simulationTimeSmallerThanSamplingInterval', ...
          'The simulation_time (%g s) is smaller than the sampling_interval (%g s).', simulation_time, sampling_interval);
end

num_samples_exact = simulation_time / sampling_interval;
num_samples_rounded = round(num_samples_exact);
if abs(num_samples_exact - num_samples_rounded) > eps
    warning('get_tppsm_data:NonIntegerRatioSamples', ...
            'simulation_time / sampling_interval is not an integer. Rounded from %.5g to %d samples.', num_samples_exact, num_samples_rounded);
end

if isempty(general_params)
    general_params = get_general_parameters();
end
general_params.simulation_time = simulation_time;
general_params.dt = sampling_interval;

if isempty(irr_params_set)
    irr_params_set = get_irregularity_parameters();
end

irr_params = irr_params_set.(scenario);

% If the external rhof_veff_ratio_L1 is not provided, compute it.
if isempty(rhof_veff_ratio_L1)
    rhof_veff_ratio_L1 = get_rhof_veff_ratio(general_params);
else
    fprintf('Using externally provided rhof_veff_ratio_L1: %g\n', rhof_veff_ratio_L1);
end

rng(seed);

[scint_field, ~, detrended_phase, ~, ~] = get_scintillation_time_series(...
    general_params, irr_params, rhof_veff_ratio_L1, seed);

psi_tppsm = scint_field(1:num_samples_rounded).';
ps_realization = detrended_phase(1:num_samples_rounded).';

end

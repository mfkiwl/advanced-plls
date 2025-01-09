function los_phase = get_los_phase(simulation_time, sampling_interval, doppler_profile)
% get_los_phase
% Generates the line-of-sight (LOS) phase time series based on a Taylor
% series expansion of any order defined by the `doppler_profile`.
%
% Syntax:
%   los_phase = get_los_phase(simulation_time, sampling_interval, doppler_profile)
%
% Description:
%   This function computes a synthetic time series for the LOS phase using
%   a Taylor series expansion of arbitrary order. The expansion coefficients
%   are provided in the `doppler_profile` vector, where each element
%   represents a coefficient for a specific order of the expansion.
%
% Inputs:
%   simulation_time    - Total duration of the simulation in seconds (positive scalar).
%   sampling_interval  - Sampling interval in seconds for the LOS phase time series (positive scalar).
%   doppler_profile    - Row vector containing Taylor series coefficients:
%                        [phase_0, fd, fdr, coeff_3, coeff_4, ..., coeff_N],
%                        where:
%                          - phase_0: Initial LOS phase (radians)
%                          - fd: Doppler frequency shift (Hz)
%                          - fdr: Doppler drift rate (Hz/s)
%                          - coeff_n: Higher-order coefficients (Hz/s^(n-1))
%
% Outputs:
%   los_phase          - Column vector containing the LOS phase time series
%                        sampled at intervals of `sampling_interval`.
%
% Notes:
%   - The LOS phase is computed using the formula:
%       los_phase(t) = phase_0 + 2*pi*sum(coeff_n * t^n/n!)
%     where:
%       - `coeff_n` is the nth-order Taylor series coefficient.
%
% Examples:
%   % Generate an LOS phase time series for a 300-second simulation, with
%   0.01 s sampling interval and a Doppler profile of [0, 1000, 0.94, 0.1].
%   doppler_profile = [0, 1000, 0.94, 0.1];
%   los_phase = get_los_phase(300, 0.01, doppler_profile);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 ORCID: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

% Input validation
validateattributes(simulation_time, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_los_phase', 'simulation_time');
validateattributes(sampling_interval, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_los_phase', 'sampling_interval');
validateattributes(doppler_profile, {'numeric'}, ...
    {'row', 'real', 'finite', 'nonnan'}, 'get_los_phase', 'doppler_profile');

if simulation_time < sampling_interval
    error('get_los_phase:simulationTimeSmallerThanSamplingInterval', ...
        'The inputed value of `simulation_time` was %g, which is smaller than the value of the `sampling_interval`, %g', ...
        simulation_time, sampling_interval)
end

% Check if simulation_time / sampling_interval is an integer
num_samples_exact = simulation_time / sampling_interval;
num_samples_rounded = round(num_samples_exact);

if abs(num_samples_exact - num_samples_rounded) > eps
    % Issue a warning if rounding was needed
    warning('get_los_phase:NonIntegerRatio', ...
            'simulation_time / sampling_interval is not an integer. The number of samples was rounded from %.5g to %d.', ...
            num_samples_exact, num_samples_rounded);
end

% Generate the time vector based on the rounded number of samples
time_vector = (1:num_samples_rounded).' * sampling_interval;

% Compute the LOS phase time series using the Taylor series expansion
% Initialize with initial phase
los_phase = doppler_profile(1) * ones(size(time_vector)); 
for n = 2:length(doppler_profile)
    los_phase = los_phase + 2 * pi * doppler_profile(n) * ...
    (time_vector.^(n-1)) / factorial(n-1);
end
end
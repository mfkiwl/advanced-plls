function los_phase = get_los_phase(simulation_time, sampling_interval, los_phase_0, fd, fdr)
% get_los_phase
% Generates the line-of-sight (LOS) phase time series based on a fixed Doppler
% frequency shift and Doppler drift using a second-order Taylor series expansion.
%
% Syntax:
%   los_phase = get_los_phase(simulation_time, sampling_interval, los_phase_0, fd, fdr)
%
% Description:
%   This function computes a synthetic time series for the LOS phase using
%   a second-order Taylor series expansion. The LOS phase is calculated for
%   each time step based on a fixed Doppler frequency shift (`fd`) and
%   Doppler drift (`fdr`), representing the dynamics of a satellite-receiver
%   geometry. The time series starts at an initial LOS phase value (`los_phase_0`)
%   and progresses over the specified simulation time.
%
% Inputs:
%   simulation_time    - Total duration of the simulation in seconds (positive scalar).
%   sampling_interval  - Sampling interval in seconds for the LOS phase time series (positive scalar).
%   los_phase_0        - Initial phase of the LOS (real scalar, in radians).
%   fd                 - Doppler frequency shift (real scalar, in Hz).
%   fdr                - Doppler drift rate (real scalar, in Hz/s).
%
% Outputs:
%   los_phase          - Column vector containing the LOS phase time series
%                        sampled at intervals of `sampling_interval`.
%
% Notes:
%   - The LOS phase is computed using the formula:
%       los_phase(t) = los_phase_0 + 2*pi*(fd * t + fdr * t^2 / 2)
%     where:
%       - `fd` is the fixed Doppler frequency shift.
%       - `fdr` is the Doppler drift rate.
%
% Examples:
%   % Generate an LOS phase time series for a 300-second simulation, with
%   0.01 s sampling interval, 1000 Hz of Doppler frequency shift and 0.94
%   Hz/s of Doppler frequency drift.
%   los_phase = get_los_phase(300, 0.01, 0, 1000, 0.94);
%
% TODO: Extend this code to enable the user to simulate the LOS dynamics
% using higher orders Taylor series expansions.
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

% Input validation
validateattributes(simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_los_phase', 'simulation_time');
validateattributes(sampling_interval, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_los_phase', 'sampling_interval');
validateattributes(los_phase_0, {'numeric'}, {'scalar', 'real', 'finite', 'nonnan'}, 'get_los_phase', 'los_phase_0');
validateattributes(fd, {'numeric'}, {'scalar', 'real', 'finite', 'nonnan'}, 'get_los_phase', 'fd');
validateattributes(fdr, {'numeric'}, {'scalar', 'real', 'finite', 'nonnan'}, 'get_los_phase', 'fdr');

if simulation_time < sampling_interval
    error('get_los_phase:simulationTimeSmallerThanSamplingInterval', ...
        'The inputed value of `simulation_time` was %g, which is smaller than the value of the `sampling_interval`, %g', simulation_time, sampling_interval)
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
time_vector = (0:num_samples_rounded-1).' * sampling_interval;

% Compute the LOS phase using vectorized operations
los_phase = los_phase_0 + 2 * pi * (fd * time_vector + fdr * (time_vector.^2) / 2);
end


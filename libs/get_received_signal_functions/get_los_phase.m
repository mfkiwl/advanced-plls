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
%   % Generate an LOS phase time series for a 10-second simulation:
%   los_phase = get_los_phase(10, 0.01, 0, 1000, 0.5);
%
% Author: Rodrigo de Lima Florindo
% Orcid: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com
% Last Modification Date: 03/01/2025

% Input validation
validate_scalar_real_positive(simulation_time, 'simulation_time');
validate_scalar_real_positive(sampling_interval, 'sampling_interval');
validate_numerical_scalar(los_phase_0, 'los_phase_0');
validate_numerical_scalar(fd, 'fd');
validate_numerical_scalar(fdr, 'fdr');

% Create the time vector
time_vector = (0:sampling_interval:simulation_time-sampling_interval).';

% Compute the LOS phase using vectorized operations
los_phase = los_phase_0 + 2 * pi * (fd * time_vector + fdr * (time_vector.^2) / 2);
end

% Helper Functions
function validate_scalar_real_positive(value, name)
    if ~isnumeric(value) || ~isscalar(value) || value <= 0
        error('get_los_phase:InvalidInput', ...
              'The input "%s" must be a positive scalar. Received: %g.', name, value);
    end
end

function validate_numerical_scalar(value, name)
    if ~isnumeric(value) || ~isscalar(value)
        error('get_los_phase:InvalidInput', ...
              'The input "%s" must be a real scalar. Received type: %s.', name, class(value));
    end
end

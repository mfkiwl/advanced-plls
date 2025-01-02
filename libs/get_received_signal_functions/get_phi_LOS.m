function phi_LOS = get_phi_LOS(simulation_time, T_I, phi_LOS_0, fd, fdr)
% get_phi_LOS
% Generates the line-of-sight (LOS) phase time series based on a fixed Doppler
% frequency shift and Doppler drift using a second-order Taylor series expansion.
%
% Syntax:
%   phi_LOS = get_phi_LOS(simulation_time, T_I, phi_LOS_0, fd, fdr)
%
% Description:
%   This function computes a synthetic time series for the LOS phase using
%   a second-order Taylor series expansion. The LOS phase is calculated for
%   each time step based on a fixed Doppler frequency shift (`fd`) and
%   Doppler drift (`fdr`), representing the dynamics of a satellite-receiver
%   geometry. The time series starts at an initial LOS phase value (`phi_LOS_0`)
%   and progresses over the specified simulation time.
%
% Inputs:
%   simulation_time - Total duration of the simulation in seconds.
%   T_I             - Sampling interval in seconds for the LOS phase time series.
%   phi_LOS_0       - Initial phase of the LOS (in radians).
%   fd              - Doppler frequency shift (Hz).
%   fdr             - Doppler drift rate (Hz/s).
%
% Outputs:
%   phi_LOS         - Column vector containing the LOS phase time series
%                     sampled at intervals of `T_I`.
%
% Notes:
%   - The LOS phase is computed using the formula:
%       phi_LOS(t) = phi_LOS_0 + 2*pi*(fd * t + fdr * t^2 / 2)
%     where:
%       - `fd` is the fixed Doppler frequency shift.
%       - `fdr` is the Doppler drift rate.
%
% Examples:
%   % Generate an LOS phase time series for a 60-second simulation with
%   % a sampling interval of 0.01 seconds, an initial LOS phase of 0 radians,
%   % a Doppler shift of 1000 Hz, and a Doppler drift rate of 0.94 Hz/s:
%   phi_LOS = get_phi_LOS(60, 0.01, 0, 1000, 0.94);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com
% Date: 01/01/2025 (Day, Month, Year)

% Preallocate the LOS phase time series vector
phi_LOS = zeros(simulation_time / T_I, 1);

% Create the time vector for the simulation
time_support = (0:T_I:simulation_time-T_I).';

% Compute the LOS phase for each time step
for k = 1:length(phi_LOS)
    phi_LOS(k) = phi_LOS_0 + 2 * pi * (fd * time_support(k) + fdr * (time_support(k)^2) / 2);
end
end
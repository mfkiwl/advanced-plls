function thermal_noise = get_thermal_noise(simulation_time, T_I, C_over_N0_dBHz)
% get_thermal_noise
% Generates additive white Gaussian noise (AWGN) to simulate thermal noise 
% in a receiver based on specified parameters.
%
% Syntax:
%   thermal_noise = get_thermal_noise(simulation_time, T_I, C_over_N0_dBHz)
%
% Description:
%   This function generates additive white Gaussian noise (AWGN) to model 
%   thermal noise in a receiver. The noise is simulated as a wide-sense 
%   stationary complex Gaussian random process. Its variance is determined 
%   based on the carrier-to-noise density ratio (C/N₀), and correlation 
%   sampling interval.
%
% Inputs:
%   simulation_time - Total duration of the simulation (seconds).
%   T_I - Integration time (seconds).
%   C_over_N0_dBHz  - Carrier-to-noise density ratio (C/N₀) in dB-Hz.
%
% Outputs:
%   thermal_noise   - Complex Gaussian noise time series with a variance 
%                     derived from the specified parameters.
%
% Notes:
%   - The continuous time noise variance is computed as:
%       σ²_η = 2 * B * (rx_mean_power / (c/n₀)) (1),
%     where c/n₀ is the linear scale equivalent of C/N₀ [dBHz] 
%     (c/n₀ = 10^(C/N₀ / 10)), and B = (IF sampling frequency) / 2.
%
%     Assuming perfect delay estimates, so that the symbol periods and the 
%     correlated samples are aligned, the equivalent complex noise variance
%     of the signal after correlation can be given by
%       σ²_η_D = σ²_η / (2 * B * T_I) (2).
%
%     The factor 2 * B * T_I represents the amount of samples that exists
%     within a correlation period.
%
%     Thus, replacing equation (1) in equation (2), and assuming without
%     loss of generality that rx_mean_power = 1, we get
%       σ²_η_D = (2 * B * (1 / (c/n₀))) / (2 * B * T_I)
%              = 1/(T_I * c/n₀).
%
%     Please, refer to [1] for further details.
%
% References
%   [1] Lopes, Rafael A. M., Felix Antreich, Friederike Fohlmeister, 
%       Martin Kriegel, and Helio K. Kuga. “Ionospheric Scintillation 
%       Mitigation With Kalman PLLs Employing Radial Basis Function 
%       Networks.” IEEE Transactions on Aerospace and Electronic Systems, 
%       2023, 1–15. https://doi.org/10.1109/TAES.2023.3281431.
%
% Example:
%   % Generate thermal noise for:
%   % - simulation_time = 600 seconds
%   % - T_I = 0.01 seconds
%   % - C/N₀ = 40 dB-Hz
%   thermal_noise = get_thermal_noise(600, 0.01, 40);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

% Input validation
validateattributes(simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'},'get_thermal_noise', 'simulation_time');
validateattributes(T_I, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'T_I');
validateattributes(C_over_N0_dBHz, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'C_over_N0_dBHz');

if simulation_time < T_I
    error('get_thermal_noise:simulationTimeSmallerThanSamplingInterval', ...
        'The inputed value of `simulation_time` was %g, which is smaller than the value of the `T_I`, %g', simulation_time, T_I)
end

% Check if simulation_time / T_I is an integer
num_samples_exact = simulation_time / T_I;
num_samples_rounded = round(num_samples_exact);

if abs(num_samples_exact - num_samples_rounded) > eps
    % Issue a warning if rounding was needed
    warning('get_thermal_noise:NonIntegerRatioSamples', ...
            'simulation_time / T_I is not an integer. The number of samples was rounded from %.5g to %d.', ...
            num_samples_exact, num_samples_rounded);
end

% Convert CN0 from dB-Hz to linear scale
c_over_n0_linear = 10^(C_over_N0_dBHz / 10);

% Compute the noise variance after integration
sigma2_eta_d = 1 / (T_I * c_over_n0_linear);

% Generate complex Gaussian noise
real_part = randn(num_samples_rounded, 1) * sqrt(sigma2_eta_d / 2);
imag_part = randn(num_samples_rounded, 1) * sqrt(sigma2_eta_d / 2);

% Combine real and imaginary parts
thermal_noise = real_part + 1j * imag_part;

end

function thermal_noise = get_thermal_noise(simulation_time, sampling_interval, rx_mean_power, C_over_N0_dBHz, B)
% get_thermal_noise
% Generates additive white Gaussian noise (AWGN) to simulate thermal noise 
% in a receiver based on specified parameters.
%
% Syntax:
%   thermal_noise = get_thermal_noise(simulation_time, sampling_interval, rx_mean_power, ...
%                                     C_over_N0_dBHz, B)
%
% Description:
%   This function generates additive white Gaussian noise (AWGN) to model 
%   thermal noise in a receiver. The noise is simulated as a wide-sense 
%   stationary complex Gaussian random process. Its variance is determined 
%   based on the carrier-to-noise density ratio (C/N₀), receiver bandwidth, 
%   and integration time.
%
% Inputs:
%   simulation_time - Total duration of the simulation (seconds).
%   sampling_interval             - Integration time (seconds).
%   rx_mean_power   - Signal power of the received signal (linear scale).
%   C_over_N0_dBHz  - Carrier-to-noise density ratio (C/N₀) in dB-Hz.
%   B               - Receiver bandwidth (Hz).
%
% Outputs:
%   thermal_noise   - Complex Gaussian noise time series with a variance 
%                     derived from the specified parameters.
%
% Notes:
%   - The noise variance is computed as:
%       σ²_η = 2 * B * (rx_mean_power / (c/n₀)),
%     where c/n₀ is the linear scale equivalent of C/N₀ [dBHz] (c/n₀ = 10^(C/N₀ / 10)).
%   - After the integrate and dump, the remaining variance after
%   normalization of the received signal using a automatic gain controller
%   (AGC) can be given by: 
%       σ²_η_discrete = σ²_η / (B * sampling_interval).
%   - The noise is generated with independent real and imaginary components.
%
% Example:
%   % Generate thermal noise for:
%   % - simulation_time = 600 seconds
%   % - sampling_interval = 0.01 seconds
%   % - rx_mean_power = 1 (linear scale)
%   % - C/N₀ = 40 dB-Hz
%   % - B = 20 MHz
%   thermal_noise = get_thermal_noise(600, 0.01, 1, 40, 2e7);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

% Input validation
validateattributes(simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'},'get_thermal_noise', 'simulation_time');
validateattributes(sampling_interval, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'sampling_interval');
validateattributes(rx_mean_power, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'rx_mean_power');
validateattributes(B, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'B');
validateattributes(C_over_N0_dBHz, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_thermal_noise', 'C_over_N0_dBHz');

if simulation_time < sampling_interval
    error('get_thermal_noise:simulationTimeSmallerThanSamplingInterval', ...
        'The inputed value of `simulation_time` was %g, which is smaller than the value of the `sampling_interval`, %g', simulation_time, sampling_interval)
end

% Check if simulation_time / sampling_interval is an integer
num_samples_exact = simulation_time / sampling_interval;
num_samples_rounded = round(num_samples_exact);

if abs(num_samples_exact - num_samples_rounded) > eps
    % Issue a warning if rounding was needed
    warning('get_thermal_noise:NonIntegerRatioSamples', ...
            'simulation_time / sampling_interval is not an integer. The number of samples was rounded from %.5g to %d.', ...
            num_samples_exact, num_samples_rounded);
end

N_int_exact = B * sampling_interval; % Number of samples per integration window
N_int_rounded = round(N_int_exact);

if abs(N_int_exact - N_int_rounded) > eps
    % Issue a warning if rounding was needed
    warning('get_thermal_noise:NonIntegerRatioAmountOfIntegrationSamples', ...
            'B * sampling_interval is not an integer. The number of samples was rounded from %.5g to %d.', ...
            N_int_exact, N_int_rounded);
end

% Convert CN0 from dB-Hz to linear scale
c_over_n0_linear = 10^(C_over_N0_dBHz / 10);

% Compute the noise variance before integration
sigma2_eta = 2 * B * (rx_mean_power / c_over_n0_linear);

% Compute the noise variance after integration
sigma2_eta_discrete = sigma2_eta / N_int_rounded;

% Generate complex Gaussian noise
real_part = randn(num_samples_rounded, 1) * sqrt(sigma2_eta_discrete / 2);
imag_part = randn(num_samples_rounded, 1) * sqrt(sigma2_eta_discrete / 2);

% Combine real and imaginary parts
thermal_noise = real_part + 1j * imag_part;

end

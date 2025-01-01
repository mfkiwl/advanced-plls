function thermal_noise = get_thermal_noise(simulation_time, T_I, rx_mean_power, C_over_N0_dBHz, B)
% get_thermal_noise
% Generates additive white Gaussian noise (AWGN) to simulate thermal noise 
% in a receiver based on specified parameters.
%
% Syntax:
%   noise = generateThermalNoise(tx_mean_power, CN0_dBHz, B, T_I, len)
%
% Description:
%   This function simulates a wide-sense stationary thermal noise in a 
%   receiver. The noise is generated as a complex Gaussian random process 
%   with variance derived from the carrier-to-noise ratio (\text{C/N}_0), 
%   receiver bandwidth, and integration time.
%
% Inputs:
%   tx_mean_power  - Signal power of the transmitted signal (linear scale).
%   C_over_N0_dBHz - Carrier-to-noise ratio in dB-Hz.
%   B              - Receiver bandwidth (Hz).
%   T_I            - Integration time (seconds).
%   len            - Length of the time series to generate.
%
% Outputs:
%   thermal_noise  - Complex white Gaussian noise time series with a
%                    specified variance.
%
% Notes:
%   - The noise variance is computed as:
%       sigma2_eta = 2 * B * (tx_mean_power / CN0)
%     and scaled for the integration time as:
%       sigma2_eta_discrete = sigma2_eta / (B * T_I).
%   - The noise is generated as a complex Gaussian process with independent
%     real and imaginary parts.
%
% Example:
%   % Generate thermal noise for a signal with:
%   % simulation_time = 600 sec, T_I = 0.01 sec, Power = 1 (linear scale), 
%   % CN0 = 40 dB-Hz, B = 20 MHz:
%   noise = get_thermal_noise(600, 0.01, 1, 40, 2e7);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com
% Date: 01/01/2025

% Convert CN0 from dB-Hz to linear scale
C_over_N0_linear = 10^(C_over_N0_dBHz / 10);

% Compute the noise variance before integration
sigma2_eta = 2 * B * (rx_mean_power / C_over_N0_linear);

% Compute the noise variance after integration
N_int = B * T_I; % Number of samples per integration window
sigma2_eta_discrete = sigma2_eta / N_int;

% Generate complex Gaussian noise
real_part = randn(simulation_time/T_I, 1) * sqrt(sigma2_eta_discrete / 2);
imag_part = randn(simulation_time/T_I, 1) * sqrt(sigma2_eta_discrete / 2);

% Combine real and imaginary parts
thermal_noise = real_part + 1j * imag_part;

end
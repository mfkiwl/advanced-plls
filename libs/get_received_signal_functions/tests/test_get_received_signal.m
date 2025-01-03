% test_get_thermal_noise
% Unit tests for the `get_thermal_noise` function, validating its
% performance in generating thermal noise with specified properties.
%
% Syntax:
%   results = runtests('test_get_thermal_noise')
%
% Description:
%   This test class ensures that the `get_thermal_noise` function meets the
%   expected behavior for generating additive white Gaussian noise (AWGN).
%   The tests validate key properties such as the mean, autocorrelation, 
%   and variance of the generated noise.
%
% Tests:
%   1. Mean Test:
%      Validates that the mean of the generated thermal noise is close to zero, 
%      ensuring no bias in the noise generation.
%
%   2. Autocorrelation Test:
%      Ensures that the autocorrelation function values of both real and
%      imaginary parts are close to 0, excluding it value at lag 0.
%
%   3. Variance Test:
%      Verifies that the variance of the generated noise matches the theoretical 
%      variance calculated from the input parameters.
%
% Inputs:
%   None directly. The test parameters (e.g., simulation time, integration 
%   time, mean power, carrier-to-noise ratio, and bandwidth) are defined within
%   the test class.
%
% Outputs:
%   None directly. The results of the tests are displayed in the MATLAB
%   Unit Test Framework output.
%
% Example:
%   % Run the test suite:
%   results = runtests('test_get_thermal_noise');
%   disp(results);
%
% Notes:
%   - This class uses the MATLAB Unit Test Framework.
%
% Author: Rodrigo de Lima Florindo
% Author's Orcid: https://orcid.org/0000-0003-0412-5583
% Author's Email: rdlfresearch@gmail.com
% Date: 03/01/2025 (Day, Month, Year)

classdef test_get_thermal_noise < matlab.unittest.TestCase
    properties
        % Test parameters
        simulation_time = 300; % Total simulation time
        T_I = 0.01; % Integration time
        rx_mean_power = 1; % Receiver mean power
        C_over_N0_dBHz = 35; % Carrier-to-noise ratio (dB-Hz)
        B = 2e7; % Bandwidth (Hz)
    end
    
    methods (TestClassSetup)
        function classSetup(testCase)
            % Add the directory containing the function to the path
            pathToAdd = fullfile(pwd, '..');
            if ~contains(path, pathToAdd)
                addpath(pathToAdd);
                testCase.addTeardown(@() rmpath(pathToAdd));
            end
        end
    end
    
    methods (Test)
        function test_mean(testCase)
            % Test if the mean of the thermal noise is close to zero
            thermal_noise = get_thermal_noise(testCase.simulation_time, testCase.T_I, ...
                                              testCase.rx_mean_power, testCase.C_over_N0_dBHz, ...
                                              testCase.B);
            testCase.assertLessThan(abs(mean(thermal_noise)), 1e-2, ...
                "The mean of the generated thermal noise is not sufficiently " + ...
                "close to zero, indicating a possible bias.");
        end
        
        function test_acf(testCase)
            % Test if the autocorrelation of thermal noise (real and imaginary parts)
            % is significant only at lag 0
            thermal_noise = get_thermal_noise(testCase.simulation_time, testCase.T_I, ...
                                              testCase.rx_mean_power, testCase.C_over_N0_dBHz, ...
                                              testCase.B);
            real_acf = autocorr(real(thermal_noise), NumLags=min(20, length(real(thermal_noise)) - 1));
            imag_acf = autocorr(imag(thermal_noise), NumLags=min(20, length(imag(thermal_noise)) - 1));
            
            testCase.verifyTrue(all(abs(real_acf(2:end)) < 5e-2), ...
                "Real part of thermal noise has significant autocorrelation beyond lag 0.");
            testCase.verifyTrue(all(abs(imag_acf(2:end)) < 5e-2), ...
                "Imaginary part of thermal noise has significant autocorrelation beyond lag 0.");
        end
        
        function test_variance(testCase)
            % Test if the variance of the generated noise matches the theoretical value
            thermal_noise = get_thermal_noise(testCase.simulation_time, testCase.T_I, ...
                                              testCase.rx_mean_power, testCase.C_over_N0_dBHz, ...
                                              testCase.B);
            variance = var(thermal_noise);

            % Compute theoretical variance
            C_over_N0_linear = 10^(testCase.C_over_N0_dBHz / 10);
            sigma2_eta = 2 * testCase.B * (testCase.rx_mean_power / C_over_N0_linear);
            N_int = testCase.B * testCase.T_I;
            theoretical_variance = sigma2_eta / N_int;

            % Check variance of generated noise
            testCase.assertLessThan(abs(variance - theoretical_variance), 1e-2, ...
                "Generated thermal noise variance differs significantly from the theoretical value.");
        end
    end
end

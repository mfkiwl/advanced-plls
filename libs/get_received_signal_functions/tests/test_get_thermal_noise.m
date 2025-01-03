classdef test_get_thermal_noise < matlab.unittest.TestCase
    % test_get_thermal_noise
    % Unit tests for the `get_thermal_noise` function, validating its ability
    % to generate additive white Gaussian noise (AWGN) accurately and handle
    % various input scenarios.
    %
    % Tests:
    %   1. Functional Tests:
    %      - Validate that the mean of the generated thermal noise is close to zero,
    %        ensuring no bias in the noise generation.
    %      - Verify that the autocorrelation function of both the real and imaginary
    %        parts is close to zero for all lags except lag 0, confirming that the
    %        noise is uncorrelated.
    %      - Ensure that the variance of the generated noise matches the theoretical
    %        variance calculated from the input parameters.
    %
    %   2. Validation Tests:
    %      - Verify the function handles edge cases and invalid inputs:
    %        * Negative or zero simulation time.
    %        * Negative or zero integration time.
    %        * Negative or zero receiver mean power.
    %        * Negative or zero bandwidth.
    %        * Non-numeric or invalid scalar inputs for the carrier-to-noise ratio.
    %
    % Example Usage:
    %   Run the test suite using MATLAB's Unit Test Framework:
    %       results = runtests('test_get_thermal_noise');
    %       disp(results);
    %
    % Author: Rodrigo de Lima Florindo
    % Orcid: https://orcid.org/0000-0003-0412-5583
    % Email: rdlfresearch@gmail.com
    % Last Modification Date: 03/01/2025

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
        %% Functional Tests
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

        %% Validation Tests
        function test_invalid_simulation_time(testCase)
            testCase.verifyError(@() get_thermal_noise(-1, testCase.T_I, testCase.rx_mean_power, ...
                                                       testCase.C_over_N0_dBHz, testCase.B), ...
                                 'get_thermal_noise:InvalidInput', ...
                                 'The function did not raise the expected error for negative simulation time.');
        end
        
        function test_invalid_integration_time(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, -0.01, testCase.rx_mean_power, ...
                                                       testCase.C_over_N0_dBHz, testCase.B), ...
                                 'get_thermal_noise:InvalidInput', ...
                                 'The function did not raise the expected error for negative integration time.');
        end
        
        function test_invalid_rx_mean_power(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, -1, ...
                                                       testCase.C_over_N0_dBHz, testCase.B), ...
                                 'get_thermal_noise:InvalidInput', ...
                                 'The function did not raise the expected error for negative mean power.');
        end
        
        function test_invalid_bandwidth(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, ...
                                                       testCase.C_over_N0_dBHz, -1), ...
                                 'get_thermal_noise:InvalidInput', ...
                                 'The function did not raise the expected error for negative bandwidth.');
        end
        
        function test_invalid_C_over_N0_dBHz(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, ...
                                                       'invalid', testCase.B), ...
                                 'get_thermal_noise:InvalidInput', ...
                                 'The function did not raise the expected error for invalid C/N0.');
        end
    end
end

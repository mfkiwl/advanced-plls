% test_get_thermal_noise
% Unit tests for the `get_thermal_noise` function, which generates additive
% white Gaussian noise (AWGN) based on specified simulation parameters.
%
% This test class ensures that the function produces noise with the correct
% statistical properties (mean, variance, and autocorrelation) and verifies
% that it appropriately handles invalid inputs.
%
% Tests:
%   - Functional Tests:
%       * Validate that the mean of the generated noise is approximately zero.
%       * Ensure that the variance matches the theoretical variance.
%       * Verify that the autocorrelation of the noise (real and imaginary
%         components) is negligible for lags other than zero.
%   - Validation Tests:
%       * Check for appropriate error handling for negative or zero values of:
%           - Simulation time.
%           - Integration time.
%           - Receiver mean power.
%           - Bandwidth.
%       * Validate behavior with non-numeric or invalid carrier-to-noise ratio.
%
% Example:
%   Run the test suite:
%       results = runtests('test_get_thermal_noise');
%       disp(results);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com
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

        %%% simulation_time

        % Non-numeric inputs
        function test_invalid_simulation_time_string(testCase)
            testCase.verifyError(@() get_thermal_noise('invalid', testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_simulation_time_logical(testCase)
            testCase.verifyError(@() get_thermal_noise(true, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_simulation_time_cell(testCase)
            testCase.verifyError(@() get_thermal_noise({42}, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        % Non-scalar inputs
        function test_invalid_simulation_time_vector(testCase)
            testCase.verifyError(@() get_thermal_noise([600, 700], testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        function test_invalid_simulation_time_matrix(testCase)
            testCase.verifyError(@() get_thermal_noise(zeros(3, 3), testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Negative values
        function test_invalid_simulation_time_negative(testCase)
            testCase.verifyError(@() get_thermal_noise(-600, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        function test_invalid_simulation_time_negative_small(testCase)
            testCase.verifyError(@() get_thermal_noise(-0.01, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Zero value
        function test_invalid_simulation_time_zero(testCase)
            testCase.verifyError(@() get_thermal_noise(0, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Non-finite values
        function test_invalid_simulation_time_infinite(testCase)
            testCase.verifyError(@() get_thermal_noise(Inf, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
            testCase.verifyError(@() get_thermal_noise(-Inf, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % NaN value
        function test_invalid_simulation_time_nan(testCase)
            testCase.verifyError(@() get_thermal_noise(NaN, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
        end

        % Empty input
        function test_invalid_simulation_time_empty(testCase)
            testCase.verifyError(@() get_thermal_noise([], testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Complex input
        function test_invalid_simulation_time_complex(testCase)
            testCase.verifyError(@() get_thermal_noise(600 + 1j, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedReal');
        end
        
        %%% T_I

        % Non-numeric inputs
        function test_invalid_T_I_string(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, 'invalid', testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_T_I_logical(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, true, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_T_I_cell(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, {0.01}, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        % Non-scalar inputs
        function test_invalid_T_I_vector(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, [0.01, 0.02], testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        function test_invalid_T_I_matrix(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, zeros(3, 3), testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Negative values
        function test_invalid_T_I_negative(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, -0.01, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Zero value
        function test_invalid_T_I_zero(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, 0, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Non-finite values
        function test_invalid_T_I_infinite(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, Inf, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, -Inf, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % NaN value
        function test_invalid_T_I_nan(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, NaN, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
        end

        % Empty input
        function test_invalid_T_I_empty(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, [], testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Complex input
        function test_invalid_T_I_complex(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, 0.01 + 1j, testCase.rx_mean_power, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedReal');
        end
        
        %%% rx_mean_power

        % Non-numeric inputs
        function test_invalid_rx_mean_power_string(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, 'invalid', testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_rx_mean_power_logical(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, true, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_rx_mean_power_cell(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, {1}, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        % Non-scalar inputs
        function test_invalid_rx_mean_power_vector(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, [1, 2], testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        function test_invalid_rx_mean_power_matrix(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, zeros(3, 3), testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Negative values
        function test_invalid_rx_mean_power_negative(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, -1, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Zero value
        function test_invalid_rx_mean_power_zero(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, 0, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Non-finite values
        function test_invalid_rx_mean_power_infinite(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, Inf, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, -Inf, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % NaN value
        function test_invalid_rx_mean_power_nan(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, NaN, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
        end

        % Empty input
        function test_invalid_rx_mean_power_empty(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, [], testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Complex input
        function test_invalid_rx_mean_power_complex(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, 1 + 1j, testCase.C_over_N0_dBHz, testCase.B), 'MATLAB:get_thermal_noise:expectedReal');
        end
        
        %%% C_over_N0_dBHz

        % Non-numeric inputs
        function test_invalid_C_over_N0_dBHz_string(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, 'invalid', testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_C_over_N0_dBHz_logical(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, true, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_C_over_N0_dBHz_cell(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, {40}, testCase.B), 'MATLAB:get_thermal_noise:invalidType');
        end

        % Non-scalar inputs
        function test_invalid_C_over_N0_dBHz_vector(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, [40, 45], testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        function test_invalid_C_over_N0_dBHz_matrix(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, zeros(3, 3), testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Non-finite values
        function test_invalid_C_over_N0_dBHz_infinite(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, Inf, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, -Inf, testCase.B), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % NaN value
        function test_invalid_C_over_N0_dBHz_nan(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, NaN, testCase.B), 'MATLAB:get_thermal_noise:expectedFinite');
        end

        % Empty input
        function test_invalid_C_over_N0_dBHz_empty(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, [], testCase.B), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Complex input
        function test_invalid_C_over_N0_dBHz_complex(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, 40 + 1j, testCase.B), 'MATLAB:get_thermal_noise:expectedReal');
        end

        %%% B

        % Non-numeric inputs
        function test_invalid_B_string(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, 'invalid'), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_B_logical(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, true), 'MATLAB:get_thermal_noise:invalidType');
        end

        function test_invalid_B_cell(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, {2e7}), 'MATLAB:get_thermal_noise:invalidType');
        end

        % Non-scalar inputs
        function test_invalid_B_vector(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, [2e7, 3e7]), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        function test_invalid_B_matrix(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, zeros(3, 3)), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Negative values
        function test_invalid_B_negative(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, -2e7), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Zero value
        function test_invalid_B_zero(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, 0), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % Non-finite values
        function test_invalid_B_infinite(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, Inf), 'MATLAB:get_thermal_noise:expectedFinite');
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, -Inf), 'MATLAB:get_thermal_noise:expectedPositive');
        end

        % NaN value
        function test_invalid_B_nan(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, NaN), 'MATLAB:get_thermal_noise:expectedFinite');
        end

        % Empty input
        function test_invalid_B_empty(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, []), 'MATLAB:get_thermal_noise:expectedScalar');
        end

        % Complex input
        function test_invalid_B_complex(testCase)
            testCase.verifyError(@() get_thermal_noise(testCase.simulation_time, testCase.T_I, testCase.rx_mean_power, testCase.C_over_N0_dBHz, 2e7 + 1j), 'MATLAB:get_thermal_noise:expectedReal');
        end
    end
end

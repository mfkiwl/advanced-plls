% test_get_los_phase
% Unit tests for the `get_los_phase` function, which computes the line-of-sight 
% (LOS) phase time series based on Doppler frequency shift and Doppler drift.
%
% This test class verifies the correctness of the LOS phase calculations 
% under various scenarios and ensures appropriate error handling for invalid inputs.
%
% Tests:
%   - Functional Tests:
%       * Validate LOS phase progression for constant Doppler shift (linear progression).
%       * Validate LOS phase progression for Doppler shift with drift (quadratic progression).
%   - Validation Tests:
%       * Check for proper error handling when:
%           - Simulation time is negative or zero.
%           - Sampling interval is negative or zero.
%           - Initial phase, Doppler shift, or Doppler drift are invalid (non-scalar or non-numeric).
%
% Example:
%   Run the test suite:
%       results = runtests('test_get_los_phase');
%       disp(results);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com
classdef test_get_los_phase < matlab.unittest.TestCase

    properties
        % Test parameters
        simulation_time = 300; % seconds
        sampling_interval = 0.01; % seconds
        los_phase_0 = 0; % radians
        fd = 1000; % Hz
        fdr = 0.94; % Hz/s
    end

    methods (TestClassSetup)
        function classSetup(testCase)
            % Add parent directory to path
            pathToAdd = fullfile(pwd, '..');
            if ~contains(path, pathToAdd)
                addpath(pathToAdd);
                testCase.addTeardown(@() rmpath(pathToAdd));
            end
        end
    end

    methods (Test)
        %% Functional Tests
        function test_constant_doppler(testCase)
            % Test LOS phase progression with constant Doppler shift (no drift)
            fd_constant = 1000; % Hz
            fdr_constant = 0; % Hz/s

            % Expected result
            time_vector = (0:testCase.sampling_interval:testCase.simulation_time-testCase.sampling_interval).';
            expected_los_phase = testCase.los_phase_0 + 2 * pi * (fd_constant * time_vector);

            % Call the function
            los_phase = get_los_phase(testCase.simulation_time, testCase.sampling_interval, ...
                                      testCase.los_phase_0, fd_constant, fdr_constant);

            % Verify results
            testCase.verifyEqual(los_phase, expected_los_phase, "AbsTol", 1e-10, ...
                "The LOS phase does not match the expected linear progression.");
        end

        function test_with_drift(testCase)
            % Test LOS phase progression with Doppler shift and drift
            time_vector = (0:testCase.sampling_interval:testCase.simulation_time-testCase.sampling_interval).';
            expected_los_phase = testCase.los_phase_0 + 2 * pi * (testCase.fd * time_vector + ...
                                      testCase.fdr * (time_vector.^2) / 2);

            % Call the function
            los_phase = get_los_phase(testCase.simulation_time, testCase.sampling_interval, ...
                                      testCase.los_phase_0, testCase.fd, testCase.fdr);

            % Verify results
            testCase.verifyEqual(los_phase, expected_los_phase, "AbsTol", 1e-10, ...
                "The LOS phase does not match the expected quadratic progression.");
        end

        %% Validation Tests
        function test_invalid_simulation_time(testCase)
            testCase.verifyError(@() get_los_phase(-1, testCase.sampling_interval, testCase.los_phase_0, testCase.fd, testCase.fdr), ...
                                 'get_los_phase:InvalidInput', ...
                                 'The function did not raise the expected error for negative simulation time.');
        end

        function test_invalid_sampling_interval(testCase)
            testCase.verifyError(@() get_los_phase(testCase.simulation_time, -0.1, testCase.los_phase_0, testCase.fd, testCase.fdr), ...
                                 'get_los_phase:InvalidInput', ...
                                 'The function did not raise the expected error for negative sampling interval.');
        end

        function test_invalid_initial_phase(testCase)
            testCase.verifyError(@() get_los_phase(testCase.simulation_time, testCase.sampling_interval, 'invalid', testCase.fd, testCase.fdr), ...
                                 'get_los_phase:InvalidInput', ...
                                 'The function did not raise the expected error for invalid initial phase.');
        end

        function test_invalid_doppler_shift(testCase)
            testCase.verifyError(@() get_los_phase(testCase.simulation_time, testCase.sampling_interval, testCase.los_phase_0, 'invalid', testCase.fdr), ...
                                 'get_los_phase:InvalidInput', ...
                                 'The function did not raise the expected error for invalid Doppler shift.');
        end

        function test_invalid_doppler_drift(testCase)
            testCase.verifyError(@() get_los_phase(testCase.simulation_time, testCase.sampling_interval, testCase.los_phase_0, testCase.fd, 'invalid'), ...
                                 'get_los_phase:InvalidInput', ...
                                 'The function did not raise the expected error for invalid Doppler drift.');
        end
    end
end

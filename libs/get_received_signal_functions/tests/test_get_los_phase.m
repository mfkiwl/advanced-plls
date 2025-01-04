% test_get_los_phase
% Unit tests for the `get_los_phase` function, which computes the line-of-sight 
% (LOS) phase time series based on Doppler frequency shift and Doppler drift.
%
% This test class verifies the correctness of the LOS phase calculations 
% under various scenarios and ensures appropriate error handling for invalid inputs.
%
% Example:
%   Run the test suite:
%       results = runtests('test_get_los_phase');
%       disp(results);
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com
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
            % Invalid simulation_time cases
            invalidCases = {
                'invalid', 'MATLAB:get_los_phase:invalidType';
                true, 'MATLAB:get_los_phase:invalidType';
                {300}, 'MATLAB:get_los_phase:invalidType';
                [300, 400], 'MATLAB:get_los_phase:expectedScalar';
                zeros(3, 3), 'MATLAB:get_los_phase:expectedScalar';
                -300, 'MATLAB:get_los_phase:expectedPositive';
                0, 'MATLAB:get_los_phase:expectedPositive';
                Inf, 'MATLAB:get_los_phase:expectedFinite';
                -Inf, 'MATLAB:get_los_phase:expectedPositive';
                NaN, 'MATLAB:get_los_phase:expectedFinite';
                [], 'MATLAB:get_los_phase:expectedScalar';
                300 + 1j, 'MATLAB:get_los_phase:expectedReal'
            };

            testCase.run_validation_tests('simulation_time', invalidCases);
        end

        function test_invalid_sampling_interval(testCase)
            % Invalid sampling_interval cases
            invalidCases = {
                'invalid', 'MATLAB:get_los_phase:invalidType';
                true, 'MATLAB:get_los_phase:invalidType';
                {0.01}, 'MATLAB:get_los_phase:invalidType';
                [0.01, 0.02], 'MATLAB:get_los_phase:expectedScalar';
                zeros(3, 3), 'MATLAB:get_los_phase:expectedScalar';
                -0.01, 'MATLAB:get_los_phase:expectedPositive';
                0, 'MATLAB:get_los_phase:expectedPositive';
                Inf, 'MATLAB:get_los_phase:expectedFinite';
                -Inf, 'MATLAB:get_los_phase:expectedPositive';
                NaN, 'MATLAB:get_los_phase:expectedFinite';
                [], 'MATLAB:get_los_phase:expectedScalar';
                0.01 + 1j, 'MATLAB:get_los_phase:expectedReal'
            };

            testCase.run_validation_tests('sampling_interval', invalidCases);
        end

        function test_invalid_los_phase_0(testCase)
            % Invalid los_phase_0 cases
            invalidCases = {
                'invalid', 'MATLAB:get_los_phase:invalidType';
                true, 'MATLAB:get_los_phase:invalidType';
                {0}, 'MATLAB:get_los_phase:invalidType';
                [0, 1], 'MATLAB:get_los_phase:expectedScalar';
                zeros(3, 3), 'MATLAB:get_los_phase:expectedScalar';
                Inf, 'MATLAB:get_los_phase:expectedFinite';
                -Inf, 'MATLAB:get_los_phase:expectedFinite';
                NaN, 'MATLAB:get_los_phase:expectedFinite';
                [], 'MATLAB:get_los_phase:expectedScalar'
            };

            testCase.run_validation_tests('los_phase_0', invalidCases);
        end

        function test_invalid_fd(testCase)
            % Invalid fd cases
            invalidCases = {
                'invalid', 'MATLAB:get_los_phase:invalidType';
                true, 'MATLAB:get_los_phase:invalidType';
                {1000}, 'MATLAB:get_los_phase:invalidType';
                [1000, 2000], 'MATLAB:get_los_phase:expectedScalar';
                zeros(2, 2), 'MATLAB:get_los_phase:expectedScalar';
                Inf, 'MATLAB:get_los_phase:expectedFinite';
                -Inf, 'MATLAB:get_los_phase:expectedFinite';
                NaN, 'MATLAB:get_los_phase:expectedFinite';
                [], 'MATLAB:get_los_phase:expectedScalar'
            };

            testCase.run_validation_tests('fd', invalidCases);
        end

        function test_invalid_fdr(testCase)
            % Invalid fdr cases
            invalidCases = {
                'invalid', 'MATLAB:get_los_phase:invalidType';
                true, 'MATLAB:get_los_phase:invalidType';
                {0.94}, 'MATLAB:get_los_phase:invalidType';
                [0.94, 1.0], 'MATLAB:get_los_phase:expectedScalar';
                zeros(3, 3), 'MATLAB:get_los_phase:expectedScalar';
                Inf, 'MATLAB:get_los_phase:expectedFinite';
                -Inf, 'MATLAB:get_los_phase:expectedFinite';
                NaN, 'MATLAB:get_los_phase:expectedFinite';
                [], 'MATLAB:get_los_phase:expectedScalar'
            };

            testCase.run_validation_tests('fdr', invalidCases);
        end
    end

    methods
        function run_validation_tests(testCase, inputName, invalidCases)
            % Iterate through invalid cases and verify errors
            for k = 1:size(invalidCases, 1)
                input = invalidCases{k, 1};
                expectedError = invalidCases{k, 2};
        
                % Generate inputs and run test
                inputs = testCase.generate_inputs(inputName, input);
                inputStr = testCase.safe_input_strings(input);
                testCase.verifyError(@() get_los_phase(inputs{:}), ...
                                     expectedError, ...
                                     sprintf('Failed for %s with input: %s', inputName, inputStr));
            end
        end
        
        function inputs = generate_inputs(testCase, fieldName, value)
            % Create inputs for the get_los_phase function
            inputs = {...
                testCase.simulation_time, ...
                testCase.sampling_interval, ...
                testCase.los_phase_0, ...
                testCase.fd, ...
                testCase.fdr};
            % Replace the appropriate field using logical indexing
            fieldNames = {'simulation_time', 'sampling_interval', 'los_phase_0', 'fd', 'fdr'};
            inputs(strcmp(fieldName, fieldNames)) = {value};
        end
        
        function str = safe_input_strings(~, input)
            % Convert input to string safely for error messages
            if ischar(input) || isstring(input)
                str = char(input);
            elseif isnumeric(input) || islogical(input)
                str = mat2str(input);
            else
                str = '<unconvertible input>';
            end
        end
    end
end
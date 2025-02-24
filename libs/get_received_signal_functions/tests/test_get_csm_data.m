classdef test_get_csm_data < matlab.unittest.TestCase
    % test_get_csm_data
    %
    % Syntax:
    %   results = runtests('test_get_csm_data')
    %
    % Description:
    %   Unit tests for get_csm_data, which generates a Cornell Scintillation 
    %   Model (CSM) time series from a struct (csm_params) containing 
    %   scintillation parameters (S4, tau0, simulation_time, sampling_interval).
    %
    %   This suite verifies:
    %   - Correct output type for valid inputs
    %   - Approximate correctness of the S4 index
    %   - Error handling for missing fields and invalid values
    %   - Warnings when simulation_time / sampling_interval is non-integer
    %
    % Example:
    %   % Run the test suite:
    %   results = runtests('test_get_csm_data');
    %   disp(results);
    %
    % Author:
    %   Rodrigo de Lima Florindo
    %   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com

    properties
        valid_csm_params % Minimal valid struct for csm_params
    end

    methods(TestClassSetup)
        function addParentPath(testCase)
            % Add parent directory and dependencies to the path
            csm_paths = genpath(fullfile(pwd,'..','..','scintillation_models/cornell_scintillation_model'));
            get_received_functions_path = fullfile(pwd,'..');
            all_paths = [csm_paths, ';' , get_received_functions_path];
            addpath(all_paths);
            testCase.addTeardown(@() rmpath(all_paths));
        end

        function setupValidParams(testCase)
            % Create a valid struct that passes all validations
            testCase.valid_csm_params = struct( ...
                'S4', 0.8, ...
                'tau0', 0.7, ...
                'simulation_time', 300, ...
                'sampling_interval', 0.01 ...
            );
        end
    end

    methods(Test)
        function testValidOutput(testCase)
            % testValidOutput
            % Ensures a valid csm_params produces a non-empty complex output.
            psi_csm = get_csm_data(testCase.valid_csm_params);
            testCase.verifyNotEmpty(psi_csm, ...
                'Output should not be empty for valid inputs.');
            testCase.verifyTrue(~isreal(psi_csm), ...
                'psi_csm should be a complex-valued array.');
        end

        function testStatisticalEvaluationOfS4(testCase)
            % testStatisticalEvaluationOfS4
            % Evaluate that the simulated S4 index is approximately correct
            % across multiple iterations.
            
            S4_target = testCase.valid_csm_params.S4;
            numIterations = 50;
            s4Differences = zeros(1, numIterations);

            for i = 1:numIterations
                psi_csm = get_csm_data(testCase.valid_csm_params);
                intensity = abs(psi_csm).^2;
                s4_calc = sqrt(var(intensity) / mean(intensity)^2); 
                % This is a quick approximate approach to the S4 index.

                s4Differences(i) = abs(s4_calc - S4_target);
            end

            meanDiff = mean(s4Differences);
            stdDiff = std(s4Differences);
            meanTolerance = 0.04; % approximate tolerance
            stdTolerance  = 0.06;

            testCase.verifyLessThanOrEqual(meanDiff, meanTolerance, ...
                sprintf('Mean difference %.4f exceeds tolerance %.4f.', meanDiff, meanTolerance));
            testCase.verifyLessThanOrEqual(stdDiff, stdTolerance, ...
                sprintf('Std difference %.4f exceeds tolerance %.4f.', stdDiff, stdTolerance));
        end

        function testS4ZeroError(testCase)
            % testS4ZeroError
            % S4 = 0 is invalid (strict positivity is required).
            invalidParams = testCase.valid_csm_params;
            invalidParams.S4 = 0;
            testCase.verifyError(@() get_csm_data(invalidParams), ...
                'MATLAB:get_csm_data:notGreater', ...
                'S4=0 should raise an expectedPositive error.');
        end

        function testTimeLessThanSamplingIntervalError(testCase)
            % testTimeLessThanSamplingIntervalError
            % If simulation_time < sampling_interval, an error is raised.
            invalidParams = testCase.valid_csm_params;
            invalidParams.simulation_time = 0.001;
            invalidParams.sampling_interval = 0.01;
            testCase.verifyError(@() get_csm_data(invalidParams), ...
                'get_csm_data:simulationTimeSmallerThanSamplingInterval', ...
                'Expected an error if simulation_time < sampling_interval.');
        end

        function testRoundingWarning(testCase)
            % testRoundingWarning
            % If simulation_time / sampling_interval is not an integer,
            % we expect a warning about rounding.
            invalid_csm_params = struct( ...
                'S4', 0.8, ...
                'tau0', 0.7, ...
                'simulation_time', 0.011, ...
                'sampling_interval', 0.01 ...
            );
            testCase.verifyWarning(@() get_csm_data(invalid_csm_params), ...
                'get_csm_data:NonIntegerRatioSamples', ...
                'Expected a rounding warning for non-integer ratio of simulation_time to sampling_interval.');
        end

        function testMissingFieldError(testCase)
            % testMissingFieldError
            % If the csm_params struct is missing one of the required fields, it raises an error.
            fields = {'S4','tau0','simulation_time','sampling_interval'};
            for f = 1:numel(fields)
                invalidParams = testCase.valid_csm_params;
                invalidParams = rmfield(invalidParams, fields{f});
                
                testCase.verifyError(@() get_csm_data(invalidParams), ...
                    'get_csm_data:MissingField', ...
                    sprintf('Expected a MissingField error when "%s" is removed.', fields{f}));
            end
        end

        function testInvalidInputType(testCase)
            % testInvalidInputType
            % If the input is not a struct, we expect a type error.
            nonStructInput = 42;
            testCase.verifyError(@() get_csm_data(nonStructInput), ...
                'MATLAB:get_csm_data:invalidType', ...
                'Expected a mustBeA error when input is not a struct.');
        end

        function testInvalidFieldValues(testCase)
            % testInvalidFieldValues
            % Check negative or nonsensical values for fields in csm_params.
            invalidParams = testCase.valid_csm_params;

            % Negative S4
            invalidParams.S4 = -0.1;
            testCase.verifyError(@() get_csm_data(invalidParams), ...
                'MATLAB:get_csm_data:notGreater', ...
                'Negative S4 should trigger an expectedPositive error.');
            
            % Reset
            invalidParams.S4 = 0.8;
            % Negative tau0
            invalidParams.tau0 = -1;
            testCase.verifyError(@() get_csm_data(invalidParams), ...
                'MATLAB:get_csm_data:expectedPositive', ...
                'Negative tau0 should trigger an expectedPositive error.');
        end
    end
end

% test_get_csm_data
% Unit tests for the `get_csm_data` function, which simulates a Cornell
% Scintillation Model (CSM) time series based on specified scintillation
% parameters.
%
% This test class ensures that the function generates a valid complex field
% time series and appropriately handles invalid input scenarios.
%
% Tests:
%   - Valid Input Tests:
%       * Ensure that valid inputs produce a non-empty output.
%   - Invalid Input Tests:
%       * Check for proper error handling when:
%           - S4 index is negative, greater than 1, non-scalar, or non-numeric.
%           - tau0 (decorrelation time) is negative, non-scalar, or non-numeric.
%           - simulation_time is negative or non-numeric.
%           - T_I (sampling time) is negative or non-numeric.
%
% Example:
%   Run the test suite:
%       results = runtests('test_get_csm_data');
%       disp(results);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com
classdef test_get_csm_data < matlab.unittest.TestCase
    
    properties
        S4 = 0.8;                   % Valid S4 index
        tau0 = 0.7;                 % Valid decorrelation time
        simulation_time = 300;      % Valid total simulation time
        T_I = 0.01;                 % Valid sampling time
    end

    methods (TestClassSetup)
        function classSetup(testCase)
            % Add parent directory to path
            pathToAdd_1 = fullfile(pwd, '..');
            if ~contains(path, pathToAdd_1)
                addpath(pathToAdd_1);
                testCase.addTeardown(@() rmpath(pathToAdd_1));
            end
            pathToAdd_2 = fullfile(pwd, '..','..','scintillation_models/cornell_scintillation_model');
            if ~contains(path, pathToAdd_2)
                addpath(pathToAdd_2);
                testCase.addTeardown(@() rmpath(pathToAdd_2));
            end
        end
    end

    methods (Test)
        %% Valid Input Tests
        function test_valid_inputs(testCase)
            % Test that valid inputs produce a result without error
            try
                psi_csm = get_csm_data(testCase.S4, testCase.tau0, ...
                                       testCase.simulation_time, testCase.T_I);
                testCase.verifyNotEmpty(psi_csm);
            catch ME
                testCase.verifyFail(['Unexpected error: ', ME.message]);
            end
        end
        
        %% Invalid Input Tests
        %%% S4 tests
        function test_invalid_S4_negative(testCase)
            % Test invalid negative S4
            testCase.verifyError(@() get_csm_data(-0.1, testCase.tau0, ...
                                                  testCase.simulation_time, testCase.T_I), ...
                                 'MATLAB:get_csm_data:notGreaterEqual');
        end
        
        function test_invalid_S4_greater_than_one(testCase)
            % Test invalid S4 greater than 1
            testCase.verifyError(@() get_csm_data(1.1, testCase.tau0, ...
                                                  testCase.simulation_time, testCase.T_I), ...
                                 'MATLAB:get_csm_data:notLessEqual');
        end

        function test_invalid_S4_non_scalar(testCase)
            % Test non-scalar S4
            testCase.verifyError(@() get_csm_data([0.8, 0.9], testCase.tau0, ...
                                                  testCase.simulation_time, testCase.T_I), ...
                                 'MATLAB:get_csm_data:expectedScalar');
        end

        function test_invalid_S4_non_numeric(testCase)
            % Test non-numeric S4
            testCase.verifyError(@() get_csm_data('invalid', testCase.tau0, ...
                                                  testCase.simulation_time, testCase.T_I), ...
                                 'MATLAB:get_csm_data:invalidType');
        end
        
        %%% tau0 tests
        function test_invalid_tau0_negative(testCase)
            % Test invalid negative tau0
            testCase.verifyError(@() get_csm_data(testCase.S4, -0.5, ...
                                                  testCase.simulation_time, testCase.T_I), ...
                                 'MATLAB:get_csm_data:expectedPositive');
        end

        function test_invalid_tau0_non_numeric(testCase)
            % Test non-numeric tau0
            testCase.verifyError(@() get_csm_data(testCase.S4, 'invalid', ...
                                                  testCase.simulation_time, testCase.T_I), ...
                                 'MATLAB:get_csm_data:invalidType');
        end

        function test_invalid_tau0_non_scalar(testCase)
            % Test non-scalar tau0
            testCase.verifyError(@() get_csm_data(testCase.S4, [0.7, 0.8], ...
                                                  testCase.simulation_time, testCase.T_I), ...
                                 'MATLAB:get_csm_data:expectedScalar');
        end
        %%% simulation_time tests
        function test_invalid_simulation_time_negative(testCase)
            % Test invalid negative simulation_time
            testCase.verifyError(@() get_csm_data(testCase.S4, testCase.tau0, ...
                                                  -10, testCase.T_I), ...
                                 'MATLAB:get_csm_data:expectedPositive');
        end
        
        function test_invalid_simulation_time_non_numeric(testCase)
            % Test non-numeric simulation_time
            testCase.verifyError(@() get_csm_data(testCase.S4, testCase.tau0, ...
                                                  'invalid', testCase.T_I), ...
                                 'MATLAB:get_csm_data:invalidType');
        end
        
        %%% T_I tests
        function test_invalid_T_I_negative(testCase)
            % Test invalid negative T_I
            testCase.verifyError(@() get_csm_data(testCase.S4, testCase.tau0, ...
                                                  testCase.simulation_time, -0.01), ...
                                 'MATLAB:get_csm_data:expectedPositive');
        end
        
        function test_invalid_T_I_non_numeric(testCase)
            % Test non-numeric T_I
            testCase.verifyError(@() get_csm_data(testCase.S4, testCase.tau0, ...
                                                  testCase.simulation_time, 'invalid'), ...
                                 'MATLAB:get_csm_data:invalidType');
        end
    end
end

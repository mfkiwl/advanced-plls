% test_get_mfpsm_data
% Unit Test Class for the `get_mfpsm_data` function.
%
% This test class verifies the functionality and robustness of the `get_mfpsm_data` 
% function, which generates Multi-Frequency Phase Screen Model (MFPSM) realizations 
% based on input scintillation parameters. The tests ensure the function produces 
% accurate results for valid inputs and handles invalid inputs gracefully.
%
% Tests:
%   - Valid Input Tests:
%       * Confirm the function produces non-empty outputs for valid inputs.
%
%   - Invalid Input Tests:
%       * Check error handling for invalid S4 values:
%           - Negative values.
%           - Values greater than 1.
%           - Non-scalar or non-numeric inputs.
%       * Verify errors for invalid tau0 inputs:
%           - Negative values.
%           - Non-scalar or non-numeric inputs.
%       * Test simulation_time for:
%           - Negative values.
%           - Non-numeric inputs.
%       * Validate T_I for:
%           - Negative values.
%           - Non-numeric inputs.
%
% Example:
%   Run the test suite:
%       results = runtests('test_get_mfpsm_data');
%       disp(results);
%
% Properties:
%   - S4: Scintillation severity index (default: 0.8).
%   - tau0: Decorrelation time in seconds (default: 0.7).
%   - simulation_time: Total simulation duration in seconds (default: 300).
%   - T_I: Sampling interval in seconds (default: 0.01).
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com
classdef test_get_mfpsm_data < matlab.unittest.TestCase
    
    properties
        S4 = 0.8;                   % Valid S4 index
        tau0 = 0.7;                 % Valid decorrelation time
        simulation_time = 300;      % Valid total simulation time
        T_I = 0.01;                 % Valid sampling time
    end

    methods (TestClassSetup)
        function classSetup(testCase)
            % Add parent directories to path
            pathToAdd_1 = fullfile(pwd, '..');
            if ~contains(path, pathToAdd_1)
                addpath(pathToAdd_1);
                testCase.addTeardown(@() rmpath(pathToAdd_1));
            end
            pathToAdd_2 = genpath(fullfile(pwd, '..', '..', ['scintillation_models' ...
                '/gnss-scintillation-simulator_2-param']));
            if ~contains(path, pathToAdd_2)
                addpath(pathToAdd_2);
                testCase.addTeardown(@() rmpath(pathToAdd_2));
            end
        end
    end

    methods (Test)
        %% Valid Input Tests
        function test_valid_inputs(testCase)
            % Test that valid inputs produce results without error
            try
                [psi_mfpsm, ps_realization] = get_mfpsm_data(testCase.S4, testCase.tau0, ...
                                                            testCase.simulation_time, testCase.T_I);
                testCase.verifyNotEmpty(psi_mfpsm);
                testCase.verifyNotEmpty(ps_realization);
            catch ME
                testCase.verifyFail(['Unexpected error: ', ME.message]);
            end
        end
        
        %% Invalid Input Tests
        %%% S4 tests
        function test_invalid_S4_negative(testCase)
            % Test invalid negative S4
            testCase.verifyError(@() get_mfpsm_data(-0.1, testCase.tau0, ...
                                                    testCase.simulation_time, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end
        
        function test_invalid_S4_greater_than_one(testCase)
            % Test invalid S4 greater than 1
            testCase.verifyError(@() get_mfpsm_data(1.1, testCase.tau0, ...
                                                    testCase.simulation_time, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end

        function test_invalid_S4_non_scalar(testCase)
            % Test non-scalar S4
            testCase.verifyError(@() get_mfpsm_data([0.8, 0.9], testCase.tau0, ...
                                                    testCase.simulation_time, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end

        function test_invalid_S4_non_numeric(testCase)
            % Test non-numeric S4
            testCase.verifyError(@() get_mfpsm_data('invalid', testCase.tau0, ...
                                                    testCase.simulation_time, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end
        
        %%% tau0 tests
        function test_invalid_tau0_negative(testCase)
            % Test invalid negative tau0
            testCase.verifyError(@() get_mfpsm_data(testCase.S4, -0.5, ...
                                                    testCase.simulation_time, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end

        function test_invalid_tau0_non_numeric(testCase)
            % Test non-numeric tau0
            testCase.verifyError(@() get_mfpsm_data(testCase.S4, 'invalid', ...
                                                    testCase.simulation_time, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end

        function test_invalid_tau0_non_scalar(testCase)
            % Test non-scalar tau0
            testCase.verifyError(@() get_mfpsm_data(testCase.S4, [0.7, 0.8], ...
                                                    testCase.simulation_time, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end
        
        %%% simulation_time tests
        function test_invalid_simulation_time_negative(testCase)
            % Test invalid negative simulation_time
            testCase.verifyError(@() get_mfpsm_data(testCase.S4, testCase.tau0, ...
                                                    -10, testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end
        
        function test_invalid_simulation_time_non_numeric(testCase)
            % Test non-numeric simulation_time
            testCase.verifyError(@() get_mfpsm_data(testCase.S4, testCase.tau0, ...
                                                    'invalid', testCase.T_I), ...
                                 'get_csm_data:InvalidInput');
        end
        
        %%% T_I tests
        function test_invalid_T_I_negative(testCase)
            % Test invalid negative T_I
            testCase.verifyError(@() get_mfpsm_data(testCase.S4, testCase.tau0, ...
                                                    testCase.simulation_time, -0.01), ...
                                 'get_csm_data:InvalidInput');
        end
        
        function test_invalid_T_I_non_numeric(testCase)
            % Test non-numeric T_I
            testCase.verifyError(@() get_mfpsm_data(testCase.S4, testCase.tau0, ...
                                                    testCase.simulation_time, 'invalid'), ...
                                 'get_csm_data:InvalidInput');
        end
    end
end
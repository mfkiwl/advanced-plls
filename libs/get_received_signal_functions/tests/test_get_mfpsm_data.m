classdef test_get_mfpsm_data < matlab.unittest.TestCase
    % test_get_mfpsm_data
    %
    % Unit Test Class for the get_mfpsm_data Function
    %
    % This class provides unit tests for the `get_mfpsm_data` function, which
    % generates Multi-Frequency Phase Screen Model (MFPSM) realizations 
    % based on scintillation parameters such as the S4 index, decorrelation
    % time (tau0), simulation duration, and sampling time.
    %
    % The tests include:
    % - Validation of the function's behavior with valid inputs.
    % - Verification of appropriate error handling for various invalid inputs.
    %
    % Example:
    %   To run the test:
    %       results = runtests('test_get_mfpsm_data');
    %
    % Properties:
    % - S4 (0.8): Valid S4 index (scintillation severity).
    % - tau0 (0.7): Valid decorrelation time in seconds.
    % - simulation_time (300): Valid simulation duration in seconds.
    % - T_I (0.01): Valid sampling time.
    %
    % Methods:
    % - test_valid_inputs: Verifies the function produces results without error
    %   for valid inputs.
    % - test_invalid_S4_*: Tests invalid S4 inputs (negative, >1, non-scalar,
    %   non-numeric).
    % - test_invalid_tau0_*: Tests invalid tau0 inputs (negative, non-scalar,
    %   non-numeric).
    % - test_invalid_simulation_time_*: Tests invalid simulation_time inputs
    %   (negative, non-numeric).
    % - test_invalid_T_I_*: Tests invalid T_I inputs (negative, non-numeric).
    %
    % Author: Rodrigo de Lima Florindo
    % Date: 03/01/2025
    % Last Modified: 03/01/2025
    % Email: rdlfresearch@gmail.com
    
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
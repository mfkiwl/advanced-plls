classdef test_get_kalman_pll_estimates < matlab.unittest.TestCase
    % test_get_kalman_pll_estimates
    %
    % This test suite verifies the behavior of get_kalman_pll_estimates.
    % It uses default values for the inputs and checks structural properties of the outputs.
    % It also tests error cases when inputs are invalid.
    
    properties
        DefaultReceivedSignal
        DefaultKalmanPLLConfig
        DefaultInitialEstimates
        DefaultTrainingModel
    end
    
    methods(TestClassSetup)
        function add_parent_path(test_case)
            % Add the necessary directory paths
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            libs_dir = fileparts(parent_dir);
            get_received_signal_functions_dir = [libs_dir,'\get_received_signal_functions'];
            tppsm_paths = genpath([libs_dir,'\scintillation_models\refactored_tppsm']);
            csm_paths = genpath([libs_dir,'\scintillation_models\cornell_scintillation_model']);
            arfit_path = [libs_dir,'\arfit'];
            addpath(parent_dir);
            addpath(get_received_signal_functions_dir);
            addpath(tppsm_paths);
            addpath(csm_paths);
            addpath(arfit_path);

            test_case.addTeardown(@() rmpath(parent_dir, get_received_signal_functions_dir, ...
                tppsm_paths, csm_paths, arfit_path));
        end

        function setDefaultValues(testCase)
            % Define a default received signal (simulate 10 time steps).
            testCase.DefaultReceivedSignal = ones(10, 1) * exp(1j*0.1);  % Complex signal
            
            % Define a default Kalman PLL config structure for training model 'CSM'
            % Using dummy matrices.
            load('dummy_kalman_pll_config.mat');
            testCase.DefaultKalmanPLLConfig = kcfg;
            
            % Define default initial estimates.
            load('dummy_initial_estimates.mat');
            testCase.DefaultInitialEstimates = initEst;
            
            % Default training model string.
            testCase.DefaultTrainingModel = 'CSM';
        end
    end

    methods(Test)
        function testValidInputs(testCase)
            % Test that valid inputs return outputs with the correct structure.
            [stateEst, errCovEst] = get_kalman_pll_estimates(testCase.DefaultReceivedSignal, ...
                testCase.DefaultKalmanPLLConfig, testCase.DefaultInitialEstimates, testCase.DefaultTrainingModel);
            
            % Check dimensions:
            % stateEstimates: [number of time steps] x [length(x_hat_init)]
            testCase.verifyEqual(size(stateEst,1), size(testCase.DefaultReceivedSignal,1));
            testCase.verifyEqual(size(stateEst,2), numel(testCase.DefaultInitialEstimates.x_hat_init));
            
            % error_covariance_estimates should be 3D with first dimension equal to number of steps.
            testCase.verifyEqual(size(errCovEst,1), size(testCase.DefaultReceivedSignal,1));
            testCase.verifyEqual(size(errCovEst,2), size(testCase.DefaultInitialEstimates.P_hat_init,1));
            testCase.verifyEqual(size(errCovEst,3), size(testCase.DefaultInitialEstimates.P_hat_init,2));
        end
        
        function testNonNumericReceivedSignal(testCase)
            % Passing a non-numeric received_signal should trigger an error.
            received_signal = 'not numeric';
            testCase.verifyError(@() get_kalman_pll_estimates(received_signal, ...
                testCase.DefaultKalmanPLLConfig, testCase.DefaultInitialEstimates, testCase.DefaultTrainingModel), ...
                'MATLAB:get_kalman_pll_estimates:invalidType');
        end
        
        function testMissingInitialEstimatesFields(testCase)
            % Remove one field from initial_estimates to trigger an error.
            badEstimates = testCase.DefaultInitialEstimates;
            badEstimates = rmfield(badEstimates, 'x_hat_init');
            testCase.verifyError(@() get_kalman_pll_estimates(testCase.DefaultReceivedSignal, ...
                testCase.DefaultKalmanPLLConfig, badEstimates, testCase.DefaultTrainingModel), ...
                'get_kalman_pll_estimates:MissingField');
        end
        
        function testInvalidTrainingModel(testCase)
            % Use a training model that does not exist in kalman_pll_config.
            invalidModel = 'TPPSM';
            testCase.verifyError(@() get_kalman_pll_estimates(testCase.DefaultReceivedSignal, ...
                testCase.DefaultKalmanPLLConfig, testCase.DefaultInitialEstimates, invalidModel), ...
                'get_kalman_pll_estimates:MissingConfigField');
        end
        
        function testInvalidTypeTrainingModel(testCase)
            % Provide a non-string for training_scint_model.
            invalidModel = 123;
            testCase.verifyError(@() get_kalman_pll_estimates(testCase.DefaultReceivedSignal, ...
                testCase.DefaultKalmanPLLConfig, testCase.DefaultInitialEstimates, invalidModel), ...
                'get_kalman_pll_estimates:InvalidType');
        end
    end
end

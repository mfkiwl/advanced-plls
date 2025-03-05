classdef test_get_kalman_pll_estimates < matlab.unittest.TestCase
    % test_get_kalman_pll_estimates
    %
    % This test suite verifies the behavior of get_kalman_pll_estimates.
    % It uses default values for the inputs and checks structural properties
    % of the outputs. It also tests error cases when inputs are invalid.
    %
    % Author:
    %   Rodrigo de Lima Florindo
    %   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com
    
    properties
        default_received_signal
        DefaultKalmanPLLConfig
        default_cache_dir
        DefaultInitialEstimates
        DefaultTrainingModel
        DefaultAdaptiveConfig
    end
    
    methods(TestClassSetup)
        function add_parent_path(test_case)
            % Add the necessary directory paths
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            libs_dir = fileparts(parent_dir);
            get_received_signal_functions_dir = fullfile(libs_dir, 'get_received_signal_functions');
            tppsm_paths = genpath(fullfile(libs_dir, 'scintillation_models', 'refactored_tppsm'));
            csm_paths = genpath(fullfile(libs_dir, 'scintillation_models', 'cornell_scintillation_model'));
            arfit_path = fullfile(libs_dir, 'arfit');
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
            testCase.default_received_signal = ones(10, 1) * exp(1j*0.1);  % Complex signal
            cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache_test_get_kalman_estimates');
            default_general_config = struct( ...
                'discrete_wiener_model_config', { {1, 3, 0.01, [0,0,1], 1} }, ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0', 0.7, 'simulation_time', 300, 'sampling_interval', 0.01, 'is_unwrapping_used', false), ...
                'C_over_N0_array_dBHz', 35, ...
                'initial_states_distributions_boundaries', { {[-pi,pi], [-5,5], [-0.1,0.1]} }, ...
                'real_doppler_profile', [0, 1000, 0.94], ...
                'augmentation_model_initializer', struct('id', 'arfit', 'model_params', struct('model_order', 3)), ...
                'is_use_cached_settings', false, ...
                'is_generate_random_initial_estimates', true ...
            );
            % Using dummy matrices
            [kcfg, initEst] = get_kalman_pll_config(default_general_config,cache_dir,false);
            % Define a default Kalman PLL config structure for training model 'CSM'
            testCase.DefaultKalmanPLLConfig = kcfg;
            % Define default initial estimates.
            testCase.DefaultInitialEstimates = initEst;
            % Default training model string.
            testCase.DefaultTrainingModel = 'CSM';
            % Define default adaptive configuration.
            testCase.DefaultAdaptiveConfig = struct('algorithm','none','hard_limited', false);
        end
    end

    methods(Test)
        function testValidInputs(testCase)
            % Test that valid inputs return outputs with the correct structure.
            [stateEst, errCovEst] = get_kalman_pll_estimates(...
                testCase.default_received_signal, ...
                testCase.DefaultKalmanPLLConfig, ...
                testCase.DefaultInitialEstimates, ...
                testCase.DefaultTrainingModel, ...
                testCase.DefaultAdaptiveConfig);
            
            % Check dimensions:
            % stateEstimates: [number of time steps] x [length(x_hat_init)]
            testCase.verifyEqual(size(stateEst,1), size(testCase.default_received_signal,1));
            testCase.verifyEqual(size(stateEst,2), numel(testCase.DefaultInitialEstimates.x_hat_init));
            
            % error_covariance_estimates should be 3D with first dimension equal to number of steps.
            testCase.verifyEqual(size(errCovEst,1), size(testCase.default_received_signal,1));
            testCase.verifyEqual(size(errCovEst,2), size(testCase.DefaultInitialEstimates.P_hat_init,1));
            testCase.verifyEqual(size(errCovEst,3), size(testCase.DefaultInitialEstimates.P_hat_init,2));
        end
        
        function testNonNumericReceivedSignal(testCase)
            % Passing a non-numeric received_signal should trigger an error.
            received_signal = 'not numeric';
            testCase.verifyError(@() get_kalman_pll_estimates(...
                received_signal, ...
                testCase.DefaultKalmanPLLConfig, ...
                testCase.DefaultInitialEstimates, ...
                testCase.DefaultTrainingModel, ...
                testCase.DefaultAdaptiveConfig), ...
                'MATLAB:get_kalman_pll_estimates:invalidType');
        end
        
        function testMissingInitialEstimatesFields(testCase)
            % Remove one field from initial_estimates to trigger an error.
            badEstimates = testCase.DefaultInitialEstimates;
            badEstimates = rmfield(badEstimates, 'x_hat_init');
            testCase.verifyError(@() get_kalman_pll_estimates(...
                testCase.default_received_signal, ...
                testCase.DefaultKalmanPLLConfig, ...
                badEstimates, ...
                testCase.DefaultTrainingModel, ...
                testCase.DefaultAdaptiveConfig), ...
                'get_kalman_pll_estimates:MissingField');
        end
        
        function testInvalidTrainingModel(testCase)
            % Use a training model that does not exist in kalman_pll_config.
            invalidModel = 'TPPSM';  % Assuming the default kcfg does not have 'TPPSM'
            testCase.verifyError(@() get_kalman_pll_estimates(...
                testCase.default_received_signal, ...
                testCase.DefaultKalmanPLLConfig, ...
                testCase.DefaultInitialEstimates, ...
                invalidModel, ...
                testCase.DefaultAdaptiveConfig), ...
                'get_kalman_pll_estimates:MissingConfigField');
        end
        
        function testInvalidTypeTrainingModel(testCase)
            % Provide a non-string for training_scint_model.
            invalidModel = 123;
            testCase.verifyError(@() get_kalman_pll_estimates(...
                testCase.default_received_signal, ...
                testCase.DefaultKalmanPLLConfig, ...
                testCase.DefaultInitialEstimates, ...
                invalidModel, ...
                testCase.DefaultAdaptiveConfig), ...
                'get_kalman_pll_estimates:InvalidType');
        end
    end
end

classdef test_get_initial_estimates < matlab.unittest.TestCase
% test_get_initial_estimates
%
% This test suite verifies the behavior of the get_initial_estimates function.
% It checks that:
%   1) When is_generate_random_initial_estimates is true, random initial estimates are generated.
%   2) When is_generate_random_initial_estimates is false, perfect estimates are used.
%   3) Required fields in general_config and kalman_pll_config are validated,
%      and errors are raised for missing or mismatched inputs.

    properties
        DefaultConfig
        DefaultKalmanPLLConfig
    end
    
    methods(TestMethodSetup)
        function setDefaults(testCase)
            % Set up a default valid configuration.
            testCase.DefaultConfig = struct( ...
                'is_generate_random_initial_estimates', true, ...
                'initial_states_distributions_boundaries', {{[-pi, pi], [-5, 5], [-0.1, 0.1]}}, ...
                'real_doppler_profile', [1000, 0.94, 0.03], ... % LOS states length = 3
                'scintillation_training_data_config', struct('scintillation_model', 'CSM') ...
            );
            
            % Set up a default kalman_pll_config for model 'CSM'
            testCase.DefaultKalmanPLLConfig = struct('CSM', struct('var_states_amount', 2, 'var_model_order', 1));
        end
        
        function addParentPath(testCase)
            parentDir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parentDir);
            testCase.addTeardown(@() rmpath(parentDir));
        end
    end
    
    methods(Test)
        function testValidRandomEstimates(testCase)
            % Test that random estimates are generated when the flag is true.
            config = testCase.DefaultConfig;
            config.is_generate_random_initial_estimates = true;
            
            initial_estimates = get_initial_estimates(config, testCase.DefaultKalmanPLLConfig);
            
            % Expected length: LOS (3) + VAR (2*1 = 2) = 5
            testCase.verifySize(initial_estimates.x_hat_init, [5, 1], 'x_hat_init size mismatch for random estimates.');
            testCase.verifySize(initial_estimates.P_hat_init, [5, 5], 'P_hat_init size mismatch for random estimates.');
        end
        
        function testValidPerfectEstimates(testCase)
            % Test that perfect estimates are generated when the flag is false.
            config = testCase.DefaultConfig;
            config.is_generate_random_initial_estimates = false;
            
            initial_estimates = get_initial_estimates(config, testCase.DefaultKalmanPLLConfig);
            
            % Expected length: LOS (3) + VAR (2*1 = 2) = 5
            testCase.verifySize(initial_estimates.x_hat_init, [5, 1], 'x_hat_init size mismatch for perfect estimates.');
            testCase.verifySize(initial_estimates.P_hat_init, [5, 5], 'P_hat_init size mismatch for perfect estimates.');
        end
        
        function testMissingConfigField(testCase)
            % Remove a required field from general_config and expect an error.
            config = testCase.DefaultConfig;
            config = rmfield(config, 'initial_states_distributions_boundaries');
            testCase.verifyError(@() get_initial_estimates(config, testCase.DefaultKalmanPLLConfig), ...
                'get_initial_estimates:MissingConfigField');
        end
        
        function testBoundaryProfileMismatch(testCase)
            % Set boundaries count different from the length of real_doppler_profile.
            config = testCase.DefaultConfig;
            config.initial_states_distributions_boundaries = {[-pi, pi], [-5, 5]}; % Only 2 boundaries but real_doppler_profile length is 3.
            testCase.verifyError(@() get_initial_estimates(config, testCase.DefaultKalmanPLLConfig), ...
                'get_initial_estimates:BoundaryProfileMismatch');
        end
        
        function testMissingKalmanSubstruct(testCase)
            % Remove the 'CSM' substruct from kalman_pll_config.
            config = testCase.DefaultConfig;
            config.is_generate_random_initial_estimates = true;
            bad_kalman = struct('TPPSM', struct(), 'none', struct());
            testCase.verifyError(@() get_initial_estimates(config, bad_kalman), ...
                'get_initial_estimates:MissingScintModelField');
        end
        
        function testMissingKalmanFields(testCase)
            % Remove a required field (e.g., var_model_order) from the 'CSM' substruct.
            config = testCase.DefaultConfig;
            bad_kalman = testCase.DefaultKalmanPLLConfig;
            bad_kalman.CSM = rmfield(bad_kalman.CSM, 'var_model_order');
            testCase.verifyError(@() get_initial_estimates(config, bad_kalman), ...
                'get_initial_estimates:MissingKalmanField');
        end
    end
end

classdef test_get_initial_estimates < matlab.unittest.TestCase
% test_get_initial_estimates
%
% This test suite verifies the behavior of the get_initial_estimates function.
% It checks that:
%   1) When is_generate_random_initial_estimates is true, random initial estimates are generated.
%   2) When is_generate_random_initial_estimates is false, perfect estimates are used.
%   3) Required fields in general_config and kalman_pll_config are validated,
%      and errors are raised for missing or mismatched inputs.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    properties
        default_general_config
        default_kalman_pll_config
    end
    
    methods(TestMethodSetup)
        function addParentPath(testCase)
            full_dir = genpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            addpath(full_dir);
            testCase.addTeardown(@() rmpath(full_dir));
        end
        function setDefaults(testCase)
            % Set up a default valid configuration.
            testCase.default_general_config = struct( ...
                'discrete_wiener_model_config', { {1, 3, 0.01, [0,0,1], 1} }, ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0', 0.7, 'simulation_time', 300, 'sampling_interval', 0.01, 'is_unwrapping_used', false), ...
                'C_over_N0_array_dBHz', 35, ...
                'initial_states_distributions_boundaries', { {[-pi,pi], [-5,5], [-0.1,0.1]} }, ...
                'real_doppler_profile', [0, 1000, 0.94], ...
                'augmentation_model_initializer', struct('id', 'arfit', 'model_params', struct('model_order', 3)), ...
                'is_use_cached_settings', false, ...
                'is_generate_random_initial_estimates', true ...
            );
            % Set up a default kalman_pll_config for model 'CSM'
            cache_file = fullfile(fileparts(mfilename('fullpath')), 'cache_test_get_initial_estimates','kalman_pll_config');
            [kalman_pll_config, is_cache_used] = get_cached_kalman_pll_config(testCase.default_general_config, cache_file, false);
            kalman_pll_config = update_cache(testCase.default_general_config, cache_file, kalman_pll_config, is_cache_used, false);
            testCase.default_kalman_pll_config = kalman_pll_config;
        end
    end
    
    methods(Test)
        function testValidRandomEstimates(testCase)
            % Test that random estimates are generated when the flag is true.
            config = testCase.default_general_config;
            config.is_generate_random_initial_estimates = true;
            
            initial_estimates = get_initial_estimates(config, testCase.default_kalman_pll_config);
            
            % Expected length: LOS (3) + VAR (2*1 = 2) = 5
            testCase.verifySize(initial_estimates.x_hat_init, [6, 1], 'x_hat_init size mismatch for random estimates.');
            testCase.verifySize(initial_estimates.P_hat_init, [6, 6], 'P_hat_init size mismatch for random estimates.');
        end
        
        function testValidPerfectEstimates(testCase)
            % Test that perfect estimates are generated when the flag is false.
            config = testCase.default_general_config;
            config.is_generate_random_initial_estimates = false;
            
            initial_estimates = get_initial_estimates(config, testCase.default_kalman_pll_config);
            
            % Expected length: LOS (3) + VAR (2*1 = 2) = 5
            testCase.verifySize(initial_estimates.x_hat_init, [6, 1], 'x_hat_init size mismatch for perfect estimates.');
            testCase.verifySize(initial_estimates.P_hat_init, [6, 6], 'P_hat_init size mismatch for perfect estimates.');
        end
        
        function testBoundaryProfileMismatch(testCase)
            % Set boundaries count different from the length of real_doppler_profile.
            config = testCase.default_general_config;
            config.initial_states_distributions_boundaries = {[-pi, pi], [-5, 5]}; % Only 2 boundaries but real_doppler_profile length is 3.
            testCase.verifyError(@() get_initial_estimates(config, testCase.default_kalman_pll_config), ...
                'get_initial_estimates:BoundaryProfileMismatch');
        end
        
        function testMissingKalmanSubstruct(testCase)
            % Remove the 'CSM' substruct from kalman_pll_config.
            config = testCase.default_general_config;
            config.is_generate_random_initial_estimates = true;
            bad_kalman = struct('TPPSM', struct(), 'none', struct());
            testCase.verifyError(@() get_initial_estimates(config, bad_kalman), ...
                'get_initial_estimates:MissingScintModelField');
        end
    end
end

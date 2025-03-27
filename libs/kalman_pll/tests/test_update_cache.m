classdef test_update_cache < matlab.unittest.TestCase
% test_update_cache
%
% This test suite verifies the behavior of the update_cache function
% after its refactoring. It checks that:
%   - With caching enabled, the function returns the input kalman_pll_config unchanged.
%   - With caching disabled, new settings are computed and saved to file.
%   - Missing or invalid inputs trigger the expected errors.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com
    
    properties
        default_general_config
        default_cache_file
        default_kalman_pll_config
        default_is_enable_cmd_print
    end
    
    methods(TestMethodSetup)
        function setupDefaults(testCase)
            testCase.default_general_config = struct( ...
                'discrete_wiener_model_config', { {1, 3, 0.01, [0,0,1], 1} }, ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0', 0.7, 'simulation_time', 300, 'sampling_interval', 0.01, 'is_unwrapping_used', false), ...
                'var_minimum_order', 1, ...
                'var_maximum_order', 6, ...
                'C_over_N0_array_dBHz', 35, ...
                'initial_states_distributions_boundaries', { {[-pi,pi], [-5,5], [-0.1,0.1]} }, ...
                'real_doppler_profile', [0, 1000, 0.94], ...
                'augmentation_model_initializer', struct('id', 'arfit', 'model_params', struct('model_order', 3)), ...
                'is_use_cached_settings', false, ...
                'is_generate_random_initial_estimates', true ...
            );
            
            % Create a temporary file name for the cache (with .mat extension)
            testCase.default_cache_file = [tempname, '.mat'];
            
            % Initialize default kalman_pll_config as an empty struct (nonempty struct)
            testCase.default_kalman_pll_config = struct();
        end
    end
    
    methods(TestClassSetup)
        function add_parent_path(test_case)
            % Add the necessary directory paths
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            libs_dir = fileparts(parent_dir);
            get_received_signal_functions_dir = fullfile(libs_dir,'get_received_signal_functions');
            tppsm_paths = genpath(fullfile(libs_dir,'scintillation_models','refactored_tppsm'));
            csm_paths = genpath(fullfile(libs_dir,'scintillation_models','cornell_scintillation_model'));
            arfit_path = fullfile(libs_dir,'arfit');
            addpath(parent_dir);
            addpath(get_received_signal_functions_dir);
            addpath(tppsm_paths);
            addpath(csm_paths);
            addpath(arfit_path);

            test_case.addTeardown(@() rmpath(parent_dir, get_received_signal_functions_dir, ...
                tppsm_paths, csm_paths, arfit_path));
        end
    end


    methods(TestMethodTeardown)
        function removeCacheFile(testCase)
            % Clean up: delete the cache file if it exists.
            if exist(testCase.default_cache_file, 'file')
                delete(testCase.default_cache_file);
            end
        end
    end
    
    methods(Test)
        function testCacheUsedReturnsSameConfig(testCase)
            % When caching is enabled, update_cache should return the input 
            % kalman_pll_config unchanged.
            is_cache_used = true;
            config = testCase.default_general_config;
            cache_file = testCase.default_cache_file;
            kalman_pll_config = testCase.default_kalman_pll_config;
            
            updated_config = update_cache(config, cache_file, kalman_pll_config, is_cache_used, testCase.default_is_enable_cmd_print);
            testCase.verifyEqual(updated_config, kalman_pll_config, ...
                'When caching is used, the output should equal the input kalman_pll_config.');
        end
        
        function testCacheNotUsedComputesNewSettings(testCase)
            % When caching is disabled, update_cache should compute new settings,
            % update kalman_pll_config, and save them to the cache file.
            is_cache_used = false;
            config = testCase.default_general_config;
            cache_file = testCase.default_cache_file;
            kalman_pll_config = testCase.default_kalman_pll_config;
            
            updated_config = update_cache(config, cache_file, kalman_pll_config, is_cache_used, testCase.default_is_enable_cmd_print);
            
            % Expect that build_kalman_pll_config adds a field named after the 
            % scintillation model (e.g., 'CSM')
            st_config = config.scintillation_training_data_config;
            modelField = st_config.scintillation_model;
            testCase.verifyTrue(isfield(updated_config, modelField), ...
                sprintf('Expected updated config to have field "%s".', modelField));
            
            % Verify that the cache file was created.
            testCase.verifyTrue(isfile(cache_file), 'Cache file should exist after computation.');
            
            % Load the saved file and compare its content to updated_config.
            data = load(cache_file, 'kalman_pll_config');
            testCase.verifyEqual(data.kalman_pll_config, updated_config, ...
                'The contents of the cache file must match the updated kalman_pll_config.');
        end
        
        function testInvalidIsCacheUsed(testCase)
            % Provide an invalid is_cache_used (non-logical) and verify that an error is thrown.
            config = testCase.default_general_config;
            cache_file = testCase.default_cache_file;
            kalman_pll_config = testCase.default_kalman_pll_config;
            is_cache_used = 'false';  % invalid type
            
            testCase.verifyError(@() update_cache(config, cache_file, kalman_pll_config, is_cache_used, testCase.default_is_enable_cmd_print), ...
                'MATLAB:update_cache:invalidType', ...
                'An error is expected when is_cache_used is not a logical scalar.');
        end
    end
end

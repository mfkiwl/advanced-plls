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
            % Define default valid inputs
            
            % LOS dynamics configuration (cell array expected by get_discrete_wiener_model)
            discrete_wiener_model_config = {1, 3, 0.01, [0, 0, 1], 1};
            
            % Scintillation training data configuration for the CSM model
            scintillation_training_data_config = struct( ...
                'scintillation_model', 'CSM', ...  % or 'TPPSM' if desired
                'S4', 0.8, ...
                'tau0', 0.7, ...
                'simulation_time', 300, ...
                'sampling_interval', 0.01 ...
            );
            
            % Other required fields for general_config
            var_minimum_order = 1;
            var_maximum_order = 6;
            C_over_N0_array_dBHz = 35;
            initial_states_distributions_boundaries = {[-pi,pi],[-5,5],[-0.05,0.05]};
            real_doppler_profile = [1000, 0.94];
            is_use_cached_settings = false;
            is_generate_random_initial_estimates = false;
            
            testCase.default_general_config = struct( ...
                'discrete_wiener_model_config', {discrete_wiener_model_config}, ...
                'scintillation_training_data_config', scintillation_training_data_config, ...
                'var_minimum_order', var_minimum_order, ...
                'var_maximum_order', var_maximum_order, ...
                'C_over_N0_array_dBHz', C_over_N0_array_dBHz, ...
                'initial_states_distributions_boundaries', {initial_states_distributions_boundaries}, ...
                'real_doppler_profile', real_doppler_profile, ...
                'is_use_cached_settings', is_use_cached_settings, ...
                'is_generate_random_initial_estimates', is_generate_random_initial_estimates ...
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
        
        function testMissingGeneralConfigField(testCase)
            % Remove a required field (e.g., discrete_wiener_model_config) from
            % general_config and verify that an error is thrown.
            config = testCase.default_general_config;
            config = rmfield(config, 'discrete_wiener_model_config');
            cache_file = testCase.default_cache_file;
            kalman_pll_config = testCase.default_kalman_pll_config;
            is_cache_used = false;
            
            testCase.verifyError(@() update_cache(config, cache_file, kalman_pll_config, is_cache_used, testCase.default_is_enable_cmd_print), ...
                'update_cache:MissingField', ...
                'An error is expected when a required field is missing from general_config.');
        end
        
        function testMissingSamplingInterval(testCase)
            % Remove the sampling_interval from scintillation_training_data_config and 
            % verify that an error is thrown.
            config = testCase.default_general_config;
            st_config = config.scintillation_training_data_config;
            st_config = rmfield(st_config, 'sampling_interval');
            config.scintillation_training_data_config = st_config;
            cache_file = testCase.default_cache_file;
            kalman_pll_config = testCase.default_kalman_pll_config;
            is_cache_used = false;
            
            testCase.verifyError(@() update_cache(config, cache_file, kalman_pll_config, is_cache_used, testCase.default_is_enable_cmd_print), ...
                'update_cache:MissingField', ...
                'An error is expected when sampling_interval is missing in scintillation_training_data_config.');
        end
        
        function testInvalidCacheFile(testCase)
            % Provide an invalid (empty) cache_file and verify that an error is thrown.
            config = testCase.default_general_config;
            cache_file = '';  % invalid cache file name
            kalman_pll_config = testCase.default_kalman_pll_config;
            is_cache_used = false;
            
            testCase.verifyError(@() update_cache(config, cache_file, kalman_pll_config, is_cache_used, testCase.default_is_enable_cmd_print), ...
                'MATLAB:update_cache:expectedNonempty', ...
                'An error is expected for an empty cache_file.');
        end
        
        function testInvalidKalmanPllConfig(testCase)
            % Provide an invalid kalman_pll_config (e.g., an empty array instead of a struct)
            config = testCase.default_general_config;
            cache_file = testCase.default_cache_file;
            kalman_pll_config = [];  % invalid kalman_pll_config
            is_cache_used = false;
            
            testCase.verifyError(@() update_cache(config, cache_file, kalman_pll_config, is_cache_used, testCase.default_is_enable_cmd_print), ...
                'MATLAB:update_cache:invalidType', ...
                'An error is expected for an invalid kalman_pll_config.');
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

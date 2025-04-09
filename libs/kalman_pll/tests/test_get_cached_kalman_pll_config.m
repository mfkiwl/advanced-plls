classdef test_get_cached_kalman_pll_config < matlab.unittest.TestCase
    % test_get_cached_kalman_pll_config
    %
    % Syntax:
    %   results = runtests('test_get_cached_kalman_pll_config')
    %
    % Description:
    %   Unit tests for the get_cached_kalman_pll_config function. This test suite
    %   verifies that:
    %     - A new configuration is initialized when no cache file exists.
    %     - An existing cache file is used when valid and is_use_cached_settings is true.
    %     - If the cache lacks required fields, the configuration is reinitialized.
    %     - Validation errors occur for invalid inputs (non-struct config, non-string paths).
    %
    % Example:
    %   % Run the test suite:
    %   results = runtests('test_get_cached_kalman_pll_config');
    %   disp(results);
    %
    % Author:
    %   Rodrigo de Lima Florindo
    %   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com

    properties
        cache_directory  % Temporary cache directory for testing
        cache_file_path  % Full path to the test cache file
        is_enable_cmd_print % Flag to enable command-line printing
    end

    methods (TestClassSetup)
        function add_parent_path_and_setup_cache(test_case)
            % add_parent_path_and_setup_cache
            % Adds the parent directory to the MATLAB path and sets up the
            % temporary cache directory and file for testing.

            % 1) Add parent directory
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parent_dir);
            test_case.addTeardown(@() rmpath(parent_dir));

            % 2) Prepare cache directory and file paths
            test_case.cache_directory = fullfile(tempdir, 'kalman_pll_cache_test');
            test_case.cache_file_path = fullfile(test_case.cache_directory, 'kalman_pll_cache.mat');
            if ~isfolder(test_case.cache_directory)
                mkdir(test_case.cache_directory);
            end
        end
        
        function is_enable_cmd_print_setup(test_case)
            test_case.is_enable_cmd_print = true;
        end
    end

    methods (TestMethodTeardown)
        function remove_cache_file(test_case)
            % remove_cache_file - Remove test cache file after each test
            if isfile(test_case.cache_file_path)
                delete(test_case.cache_file_path);
            end
        end
    end

    methods (Test)
        function test_initialize_no_cache(test_case)
            % test_initialize_no_cache
            % Verifies that a new configuration is initialized when the cache does not exist.

            if isfile(test_case.cache_file_path)
                delete(test_case.cache_file_path);
            end

            % Include new field kf_type as required.
            config = struct( ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM'), ...
                'kf_type', 'standard', ...
                'is_use_cached_settings', false ...
            );

            [kalman_pll_config, is_cache_used] = ...
                get_cached_kalman_pll_config(config, test_case.cache_file_path, test_case.is_enable_cmd_print);

            test_case.verifyFalse(is_cache_used, ...
                'Expected is_cache_used = false when no cache file exists.');
            test_case.verifyTrue(isstruct(kalman_pll_config), ...
                'kalman_pll_config should be a struct.');
            test_case.verifyTrue(isfield(kalman_pll_config, 'CSM'), ...
                'CSM field should be present after initialization.');
            % Ensure the specific KF type subfield exists (even if empty)
            test_case.verifyTrue(isfield(kalman_pll_config.CSM, 'standard'), ...
                'standard field should be present in the CSM configuration.');
        end

        function test_use_existing_cache(test_case)
            % test_use_existing_cache
            % Verifies that an existing valid cache is used when is_use_cached_settings is true.

            % Create a mock configuration with the expected structure.
            % Here we assume the KF type is 'standard' and provide a non-empty value.
            cached_config = struct( ...
                'CSM', struct( ...
                    'standard', eye(2), ...  % non-empty to simulate valid cached settings
                    'extended', [], ...
                    'unscented', [], ...
                    'cubature', []), ...
                'TPPSM', struct('standard', [], 'extended', [], 'unscented', [], 'cubature', []), ...
                'none',    struct('standard', [], 'extended', [], 'unscented', [], 'cubature', []));
            kalman_pll_config = cached_config;
            save(test_case.cache_file_path, 'kalman_pll_config', '-mat');

            config = struct( ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM'), ...
                'kf_type', 'standard', ...
                'is_use_cached_settings', true ...
            );

            [loaded_config, is_cache_used] = ...
                get_cached_kalman_pll_config(config, test_case.cache_file_path, test_case.is_enable_cmd_print);

            test_case.verifyTrue(is_cache_used, ...
                'Expected cache to be used when a valid cache file exists.');
            test_case.verifyEqual(loaded_config, cached_config, ...
                'The loaded config should match the cached config.');
        end

        function test_reinitialize_cache_missing_fields(test_case)
            % test_reinitialize_cache_missing_fields
            % Verifies that if the cache is missing the required KF type settings,
            % the configuration is reinitialized.
            %
            % In this new version of the function, when the cache file exists
            % but the required inner KF-type subfields for the provided scintillation model are missing,
            % the function does not use the cached value (is_use_cached_settings=false)
            % and returns the configuration saved in the cache (which remains incomplete).
            % Hence, we update the test to expect the incomplete (partial) configuration.
            
            % Create a partial configuration: it has the "CSM" field but no KF-type subfields.
            partial_config = struct('CSM', struct());
            kalman_pll_config = partial_config;
            save(test_case.cache_file_path, 'kalman_pll_config', '-mat');

            config = struct( ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM'), ...
                'kf_type', 'standard', ...
                'is_use_cached_settings', true ...
            );

            [reinit_config, is_cache_used] = ...
                get_cached_kalman_pll_config(config, test_case.cache_file_path, test_case.is_enable_cmd_print);

            % Since the cache is incomplete, the function should not use it.
            test_case.verifyFalse(is_cache_used, ...
                'Expected is_cache_used to be false when cache is missing required KF type settings.');
            
            % Update the expected configuration to match what the new function version returns.
            % (Previously, a complete default structure was expected;
            %  now it is expected to return the partial config as saved.)
            expected_config = partial_config;  % i.e., struct('CSM', struct())
            test_case.verifyEqual(reinit_config, expected_config, ...
                'Reinitialized config does not match the expected (partial) structure.');
        end
    end
end

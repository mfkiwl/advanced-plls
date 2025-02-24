classdef test_get_cached_kalman_pll_config < matlab.unittest.TestCase
% get_cached_kalman_pll_config
%
% Syntax:
%   results = runtests('get_cached_kalman_pll_config')
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

            config = struct( ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM'), ...
                'is_use_cached_settings', false ...
            );

            [kalman_pll_config, is_cache_used] = ...
                get_cached_kalman_pll_config(config, test_case.cache_file_path);

            test_case.verifyFalse(is_cache_used, ...
                'Expected is_cache_used = false when no cache file exists.');
            test_case.verifyTrue(isstruct(kalman_pll_config), ...
                'kalman_pll_config should be a struct.');
            test_case.verifyTrue(isfield(kalman_pll_config, 'CSM'), ...
                'CSM field should be present after initialization.');
        end

        function test_use_existing_cache(test_case)
            % test_use_existing_cache
            % Verifies that an existing valid cache is used when is_use_cached_settings = true.

            % Create a mock config in the cache
            cached_config = struct('CSM', struct('F', eye(2)));
            kalman_pll_config = cached_config;
            save(test_case.cache_file_path, 'kalman_pll_config', '-mat');

            config = struct( ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM'), ...
                'is_use_cached_settings', true ...
            );

            [loaded_config, is_cache_used] = ...
                get_cached_kalman_pll_config(config, test_case.cache_file_path);

            test_case.verifyTrue(is_cache_used, ...
                'Expected cache to be used when a valid cache file exists.');
            test_case.verifyEqual(loaded_config, cached_config, ...
                'The loaded config should match the cached config.');
        end

        function test_reinitialize_cache_missing_fields(test_case)
            % test_reinitialize_cache_missing_fields
            % Verifies that if the cache is missing required fields,
            % the configuration is reinitialized.

            partial_config = struct('dummy_field', 123);
            kalman_pll_config = partial_config;
            save(test_case.cache_file_path, 'kalman_pll_config', '-mat');

            config = struct( ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM'), ...
                'is_use_cached_settings', true ...
            );

            [reinit_config, is_cache_used] = ...
                get_cached_kalman_pll_config(config, test_case.cache_file_path);

            test_case.verifyFalse(is_cache_used, ...
                'Expected a new config if the cache lacks the requested model field.');
            test_case.verifyTrue(isfield(reinit_config, 'CSM'), ...
                'CSM field should be present after re-initialization.');
        end

        function test_invalid_config_type(test_case)
            % test_invalid_config_type
            % Verifies that a non-struct config triggers an error.

            invalid_config = 123; % Not a struct
            test_case.verifyError(...
                @() get_cached_kalman_pll_config(invalid_config, test_case.cache_file_path), ...
                'MATLAB:get_cached_kalman_pll_config:invalidType', ...
                'Expected an error for invalid config type.');
        end

        function test_invalid_cache_path(test_case)
            % test_invalid_cache_path
            % Verifies that a non-string path triggers an error.

            config = struct( ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM'), ...
                'is_use_cached_settings', false ...
            );
            invalid_path = 42; % Not a string or char

            test_case.verifyError(...
                @() get_cached_kalman_pll_config(config, invalid_path), ...
                'MATLAB:fileparts:MustBeChar', ...
                'Expected an error for invalid cache path type.');
        end
    end
end

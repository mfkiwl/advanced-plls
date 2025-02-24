classdef test_build_kalman_pll_config < matlab.unittest.TestCase
    % test_build_kalman_pll_config
    %
    % Syntax:
    %   results = runtests('test_build_kalman_pll_config')
    %
    % Description:
    %   Unit tests for the build_kalman_pll_config function. This test suite verifies:
    %   - Proper validation of required config fields.
    %   - Behavior when is_cache_used is true or false.
    %   - Handling of 'none' vs. 'CSM' scintillation models.
    %   - Integration with get_discrete_wiener_model and compute_settings calls.
    %
    % Example:
    %   % Run the test suite:
    %   results = runtests('test_build_kalman_pll_config');
    %   disp(results);
    %
    % Author:
    %   Rodrigo de Lima Florindo
    %   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com

    properties
        default_cache_file = 'kalman_pll_cache.mat';  % Default path used by the function
        default_sampling_interval = 0.01;
        minimal_pll_config        % Minimal struct for kalman_pll_config input
        valid_config              % Minimal valid config for multiple tests
    end

    methods(TestClassSetup)
        function add_parent_path(test_case)
            % Add the parent directory so that build_kalman_pll_config.m 
            % and its dependencies can be found.
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parent_dir);
            test_case.addTeardown(@() rmpath(parent_dir));
        end

        function setup_test_data(test_case)
            % Setup minimal structures for testing
            
            % Minimal existing kalman_pll_config
            test_case.minimal_pll_config = struct();
            
            % A valid config, referencing only mandatory fields
            % Additional fields can be appended in test methods if needed.
            test_case.valid_config = struct( ...
                'training_scint_model', 'CSM', ...
                'discrete_wiener_model_config', { {1,3,0.01,[0,0,1],1} }, ...
                'scintillation_training_data_config', { {0.8,0.7,300,0.01} }, ...
                'var_minimum_order', 1, ...
                'var_maximum_order', 6, ...
                'C_over_N0_array_dBHz', 35, ...
                'is_refractive_effects_removed', false ...
            );
        end
    end

    methods(Test)
        function test_valid_config_csm(test_case)
            % test_valid_config_csm
            % Uses a 'CSM' scintillation model and is_cache_used = false.
            config = test_case.valid_config;
            is_cache_used = false;

            updated_pll_config = build_kalman_pll_config(config, ...
                test_case.default_cache_file, test_case.default_sampling_interval, ...
                test_case.minimal_pll_config, is_cache_used);

            test_case.verifyTrue(isstruct(updated_pll_config), ...
                'Expected updated_pll_config to be a struct.');
            test_case.verifyTrue(isfield(updated_pll_config, 'CSM'), ...
                'Expected a "CSM" field in the returned config.');
        end

        function test_valid_config_none(test_case)
            % test_valid_config_none
            % Uses 'none' model to ensure only LOS states are used.
            config = test_case.valid_config;
            config.training_scint_model = 'none';
            is_cache_used = false;

            updated_pll_config = build_kalman_pll_config(config, ...
                test_case.default_cache_file, test_case.default_sampling_interval, ...
                test_case.minimal_pll_config, is_cache_used);

            test_case.verifyTrue(isstruct(updated_pll_config), ...
                'Expected a struct for updated_pll_config.');
            test_case.verifyTrue(isfield(updated_pll_config, 'none'), ...
                'Expected a "none" field in the returned config.');
            test_case.verifyEmpty(updated_pll_config.none.F_var, ...
                'F_var should be empty for "none" model.');
            test_case.verifyEmpty(updated_pll_config.none.Q_var, ...
                'Q_var should be empty for "none" model.');
        end

        function test_missing_required_field(test_case)
            % test_missing_required_field
            % Remove a required config field (e.g., is_refractive_effects_removed).
            config = rmfield(test_case.valid_config, 'is_refractive_effects_removed');
            is_cache_used = false;

            test_case.verifyError(@() build_kalman_pll_config(config, ...
                test_case.default_cache_file, test_case.default_sampling_interval, ...
                test_case.minimal_pll_config, is_cache_used), ...
                'build_kalman_pll_config:MissingField', ...
                'Expected an error due to a missing required config field.');
        end

        function test_invalid_config_type(test_case)
            % test_invalid_config_type
            % Provide a non-struct config to check validation.
            invalid_config = 123;
            is_cache_used = false;

            test_case.verifyError(@() build_kalman_pll_config(invalid_config, ...
                test_case.default_cache_file, test_case.default_sampling_interval, ...
                test_case.minimal_pll_config, is_cache_used), ...
                'MATLAB:validators:mustBeA', ...
                'Expected an error for invalid config type.');
        end

        function test_cache_used_flag(test_case)
            % test_cache_used_flag
            % If is_cache_used is true, it should skip recomputing.
            % We won't do advanced mocking of the actual caching logic, 
            % just ensure no error is thrown and the function runs.
            config = test_case.valid_config;
            is_cache_used = true;

            updated_pll_config = build_kalman_pll_config(config, ...
                test_case.default_cache_file, test_case.default_sampling_interval, ...
                test_case.minimal_pll_config, is_cache_used);

            test_case.verifyTrue(isstruct(updated_pll_config), ...
                'Should still return a struct when is_cache_used=true.');
            % We can also check the console output or use a mock approach, 
            % but for now we just confirm it didn't blow up.
        end
    end
end

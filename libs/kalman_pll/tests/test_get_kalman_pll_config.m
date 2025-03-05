classdef test_get_kalman_pll_config < matlab.unittest.TestCase
    % test_get_kalman_pll_config
    %
    % This test suite verifies the behavior of the get_kalman_pll_config function.
    % It checks that:
    %   - With a valid configuration, the function returns a configuration struct
    %     with the expected model field (added by the dummy update_cache) and a
    %     column vector of initial estimates.
    %   - Missing required fields, mismatched sampling intervals, invalid boundaries,
    %     an empty real_doppler_profile, or negative C/N0 values trigger errors.
    %
    % Dummy lower-level functions (get_cached_kalman_pll_config, update_cache,
    % and get_initial_estimates) are overridden by static methods defined in this class.
    %
    % Author:
    %   Rodrigo de Lima Florindo
    %   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com

    properties
        default_config
        default_cache_dir
        default_is_enable_cmd_print
    end

    methods(TestClassSetup)
        function add_parent_path(test_case)
            % Add the necessary directory paths.
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
            test_case.addTeardown(@() rmpath(parent_dir, get_received_signal_functions_dir, tppsm_paths, csm_paths, arfit_path));
        end
    end

    methods(TestMethodSetup)
        function setDefaultConfig(testCase)
            % Define a default valid configuration struct.
            testCase.default_config = struct( ...
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
            testCase.default_cache_dir = fullfile(fileparts(mfilename('fullpath')),'cache');
            testCase.default_is_enable_cmd_print = true;
        end
    end

    methods(Test)
        function testValidConfig(testCase)
            % Using the default configuration, verify that the outputs are valid.
            [kcfg, initEst] = get_kalman_pll_config(testCase.default_config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print);
            testCase.verifyTrue(isstruct(kcfg), 'kalman_pll_config must be a struct.');
            % Dummy update_cache adds field "CSM" to simulate a valid configuration.
            testCase.verifyTrue(isfield(kcfg, 'CSM'), 'Expected field "CSM" not found.');
            
            testCase.verifyTrue(isstruct(initEst), 'initial_estimates must be numeric.');
            testCase.verifyEqual(size(initEst,2), 1, 'initial_estimates must be a column vector.');
        end
        
        function testMissingRequiredField(testCase)
            % Remove a required field from the scintillation_training_data_config.
            config = testCase.default_config;
            config.scintillation_training_data_config = rmfield(config.scintillation_training_data_config, 'scintillation_model');
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:nonExistentField');
        end
        
        function testSamplingIntervalMismatch(testCase)
            % Modify the sampling interval in discrete_wiener_model_config so that it mismatches.
            config = testCase.default_config;
            config.discrete_wiener_model_config{3} = 0.02; % Expected is 0.01 as in scintillation_training_data_config.
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'get_kalman_pll_config:SamplingIntervalMismatch');
        end
        
        function testInvalidBoundaries(testCase)
            % Provide an invalid boundary (min >= max) for one element.
            config = testCase.default_config;
            config.initial_states_distributions_boundaries{3} = [0.1, -0.1];
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'get_kalman_pll_config:InvalidBoundaries');
        end
        
        function testEmptyRealDopplerProfile(testCase)
            % Set real_doppler_profile to an empty array.
            config = testCase.default_config;
            config.real_doppler_profile = [];
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:get_kalman_pll_config:expectedNonempty');
        end
        
        function testNegativeCOverN0(testCase)
            % Set a negative value for C_over_N0_array_dBHz.
            config = testCase.default_config;
            config.C_over_N0_array_dBHz = -35;
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:get_kalman_pll_config:expectedPositive');
        end

        function test_augmentation_model_initializer_validation(testCase)
            % Test invalid augmentation_model_initializer inputs.
            config = testCase.default_config;
            
            % Case 1: Empty initializer.
            config.augmentation_model_initializer.id = [];
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:undefinedVarOrClass');
            
            % Case 2: Non-string initializer.
            config.augmentation_model_initializer.id = 3;
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:undefinedVarOrClass'); % Error id from validateattributes.
            
            % Case 3: String not among allowed values.
            config.augmentation_model_initializer.id = 'invalid';
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:undefinedVarOrClass');
        end
    
        function testAugmentationModelMissingModelParams(testCase)
            % Remove the model_params field altogether.
            config = testCase.default_config;
            config.augmentation_model_initializer = rmfield(config.augmentation_model_initializer, 'model_params');
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:nonExistentField');
        end
    
        function testAugmentationModelEmptyModelParams(testCase)
            % Set model_params to an empty struct.
            config = testCase.default_config;
            config.augmentation_model_initializer.model_params = struct();
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:nonExistentField');
        end
    
        function testAugmentationModelMissingSubfieldForArfit(testCase)
            % For 'arfit', remove the required field 'model_order'.
            config = testCase.default_config;
            config.augmentation_model_initializer.id = 'arfit';
            config.augmentation_model_initializer.model_params = rmfield(config.augmentation_model_initializer.model_params, 'model_order');
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:nonExistentField');
        end
    
        function testAugmentationModelMissingSubfieldForAryule(testCase)
            % For 'aryule', remove the required field 'model_order'.
            config = testCase.default_config;
            config.augmentation_model_initializer.id = 'aryule';
            config.augmentation_model_initializer.model_params = rmfield(config.augmentation_model_initializer.model_params, 'model_order');
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:nonExistentField');
        end
    
        function testAugmentationModelMissingSubfieldForRbf(testCase)
            % For 'rbf', remove the required field 'neurons_amount'.
            config = testCase.default_config;
            config.augmentation_model_initializer.id = 'rbf';
            % Provide a model_params struct without 'neurons_amount'.
            config.augmentation_model_initializer.model_params = struct();
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'MATLAB:RBFUnavailable');
        end
    
        function testValidAugmentationModelNone(testCase)
            % For initializer 'none', as long as model_params is a nonempty struct, the function should succeed.
            config = testCase.default_config;
            config.augmentation_model_initializer.id = 'none';
            config.augmentation_model_initializer.model_params = struct('dummy', 1); % Dummy nonempty struct.
            [kcfg, initEst] = get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print);
            testCase.verifyTrue(isstruct(kcfg));
            testCase.verifyTrue(isstruct(initEst));
        end
    
        function testValidAugmentationModelAryule(testCase)
            % For a valid 'aryule' initializer.
            config = testCase.default_config;
            config.augmentation_model_initializer.id = 'aryule';
            config.augmentation_model_initializer.model_params = struct('model_order', 4);
            [kcfg, ~] = get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print);
            testCase.verifyTrue(isstruct(kcfg));
        end
        function testInvalidAugmentationModelAryuleWithMultiFrequency(testCase)
            % For a valid 'aryule' initializer.
            config = testCase.default_config;
            config.discrete_wiener_model_config{1} = 2;
            config.discrete_wiener_model_config{4} = [0,0,0,1];
            config.discrete_wiener_model_config{5} = [1,0.9];
            config.augmentation_model_initializer.id = 'aryule';
            config.augmentation_model_initializer.model_params = struct('model_order', 4);
            testCase.verifyError(@() get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print), ...
                'get_kalman_pll_config:incompatible_model_with_multi_frequency_tracking');
        end
        
        % This test should be uncommented when RBF case is being developed.
        % function testValidAugmentationModelRbf(testCase)
        %     % For a valid 'rbf' initializer.
        %     config = testCase.default_config;
        %     config.augmentation_model_initializer.id = 'rbf';
        %     config.augmentation_model_initializer.model_params = struct('neurons_amount', 10);
        %     [kcfg, ~] = get_kalman_pll_config(config, testCase.default_cache_dir, testCase.default_is_enable_cmd_print);
        %     testCase.verifyTrue(isstruct(kcfg));
        % end
    end
    
    methods(Static)
        % Dummy implementations to override lower-level functions for testing.
        function [kcfg, is_cache_used] = get_cached_kalman_pll_config(~, ~, ~)
            % Return an empty config and false flag to simulate no cached configuration.
            kcfg = struct();
            is_cache_used = false;
        end
        
        function kcfg = update_cache(~, ~, kcfg, ~, ~)
            % Return a dummy configuration struct with a field "CSM" to simulate a valid update.
            kcfg.CSM = struct( ...
                'F', eye(2), 'Q', eye(2), 'H', ones(1,2), 'R', eye(2), 'W', [0,0], 'augmentation_model_initializer', struct());
        end
        
        function estimates = get_initial_estimates(~, kcfg)
            % Return a dummy column vector for initial estimates.
            if isfield(kcfg, 'CSM')
                n = size(kcfg.CSM.F, 1);
            else
                n = 2;
            end
            estimates = ones(n, 1);
        end
    end
end

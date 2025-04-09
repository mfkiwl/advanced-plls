classdef test_build_kalman_pll_config < matlab.unittest.TestCase
% test_build_kalman_pll_config
%
% This test suite verifies the behavior of the build_kalman_pll_config
% function. It checks that valid inputs execute successfully and that
% invalid inputs raise appropriate errors or warnings.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    properties
        default_kalman_pll_config
        default_general_config
        default_F_los
        default_Q_los
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

    methods(TestMethodSetup)
        function setupDefaults(testCase)
            % Set up default valid inputs
            testCase.default_kalman_pll_config = struct();
            testCase.default_general_config = struct( ...
                'kf_type', 'standard', ...
                'discrete_wiener_model_config', { {1, 3, 0.01, [0,0,1], 1} }, ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0', 0.7, 'simulation_time', 300, 'sampling_interval', 0.01, 'is_unwrapping_used', false), ...
                'C_over_N0_array_dBHz', 35, ...
                'initial_states_distributions_boundaries', { {[-pi,pi], [-5,5], [-0.1,0.1]} }, ...
                'real_doppler_profile', [0, 1000, 0.94], ...
                'augmentation_model_initializer', struct('id', 'arfit', 'model_params', struct('model_order', 3)), ...
                'is_use_cached_settings', false, ...
                'is_generate_random_initial_estimates', true ...
            );
            testCase.default_F_los = eye(3);
            testCase.default_Q_los = eye(3);
        end
    end

    methods(Test)
        function testValidInputs(testCase)
            % Verify that valid inputs (default 'arfit') execute without error
            kalman_pll_config = build_kalman_pll_config( ...
                testCase.default_general_config, ...
                testCase.default_kalman_pll_config);
            
            % Check that the output struct has a field named after the scintillation model.
            modelName = testCase.default_general_config.scintillation_training_data_config.scintillation_model;
            testCase.verifyTrue(isfield(kalman_pll_config, modelName), ...
                sprintf('Expected the output struct to have a field "%s".', modelName));
            
            % Check that the substruct has all the required fields.
            requiredFields = {'F','Q','H','R','W', 'augmentation_model_initializer'};
            configSubstruct = kalman_pll_config.(modelName);
            for i = 1:length(requiredFields)
                testCase.verifyTrue(isfield(configSubstruct, requiredFields{i}), ...
                    sprintf('Missing required field "%s" in output struct.', requiredFields{i}));
            end
        end
        
        function testValidInputsAryule(testCase)
            % Verify that valid inputs using the 'aryule' initializer execute without error.
            config = testCase.default_general_config;
            config.augmentation_model_initializer.id = 'aryule';
            config.augmentation_model_initializer.model_params = struct('model_order', 4);
            
            kalman_pll_config = build_kalman_pll_config(config, testCase.default_kalman_pll_config);
            modelName = config.scintillation_training_data_config.scintillation_model;
            testCase.verifyTrue(isfield(kalman_pll_config, modelName), ...
                sprintf('Expected the output struct to have a field "%s".', modelName));
        end
        
        function testInvalidAugmentationModel(testCase)
            % Verify that an invalid augmentation_model_initializer id leads to an error.
            config = testCase.default_general_config;
            config.augmentation_model_initializer.id = 'invalid';
            % With an invalid id, none of the switch cases match and required variables remain undefined.
            testCase.verifyError(@() build_kalman_pll_config(config, testCase.default_kalman_pll_config), ...
                'MATLAB:UndefinedAugmentationModel', 'Expected error for invalid augmentation_model_initializer id was not thrown.');
        end
        
        function testDifferentScintillationModel(testCase)
            % Verify that using a different scintillation model name correctly stores results
            config = testCase.default_general_config;
            config.scintillation_training_data_config.scintillation_model = 'CSM';
            % Using a valid initializer branch (aryule in this example).
            config.augmentation_model_initializer.id = 'aryule';
            config.augmentation_model_initializer.model_params = struct('model_order', 4);
            
            kalman_pll_config = build_kalman_pll_config(config, testCase.default_kalman_pll_config);
            testCase.verifyTrue(isfield(kalman_pll_config, 'CSM'), ...
                'Expected the output struct to have a field "CSM".');
        end
    end
end
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

    methods(TestMethodSetup)
        function setupDefaults(testCase)
            % Set up default valid inputs
            testCase.default_kalman_pll_config = struct();
            testCase.default_general_config = struct( ...
                'discrete_wiener_model_config', { {1, 3, 0.01, [0,0,1], 1} }, ...
                'scintillation_training_data_config', struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0', 0.7, 'simulation_time', 300, 'sampling_interval', 0.01), ...
                'var_minimum_order', 1, ...
                'var_maximum_order', 6, ...
                'C_over_N0_array_dBHz', 35, ...
                'initial_states_distributions_boundaries', { {[-pi,pi], [-5,5], [-0.1,0.1]} }, ...
                'real_doppler_profile', [0, 1000, 0.94], ...
                'augmentation_model_initializer', 'arfit', ...
                'is_use_cached_settings', false, ...
                'is_generate_random_initial_estimates', true ...
            );
            testCase.default_F_los = eye(3);
            testCase.default_Q_los = eye(3);
        end
    end

    methods(Test)
        function testValidInputs(testCase)
            % Verify that valid inputs execute without error
            kalman_pll_config = build_kalman_pll_config( ...
                testCase.default_general_config, ...
                testCase.default_kalman_pll_config);
            
            % Check that the function adds a field named after the model
            testCase.verifyTrue(isfield(kalman_pll_config, ...
                testCase.default_general_config.scintillation_training_data_config.scintillation_model), ...
                'Expected the output struct to have a field for the model.');
        end
    end
end

classdef test_get_augmentation_model < matlab.unittest.TestCase
    properties
        training_data
        general_config
        sampling_interval
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
        function setup(testCase)
            % Define some dummy training data.
            testCase.training_data = (1:10)'; 
            testCase.sampling_interval = 0.01;
            % Set up a basic general_config.
            testCase.general_config.C_over_N0_array_dBHz = 35;
        end
    end
    
    methods(Test)
        function testArfitBranch(testCase)
            % --- Test ARFIT branch ---
            testCase.general_config.augmentation_model_initializer.id = 'arfit';
            testCase.general_config.augmentation_model_initializer.model_params.model_order = 2;
            
            % The local stubs for arfit and construct_var_matrices (defined below)
            % will provide predictable outputs.
            augData = get_augmentation_model(testCase.training_data, testCase.general_config, testCase.sampling_interval);
            
            testCase.verifyEqual(lower(augData.augmentation_type), 'arfit', 'Augmentation type should be arfit.');
            testCase.verifyTrue(isfield(augData, 'F_aug'), 'F_aug field must exist.');
            testCase.verifyTrue(isfield(augData, 'Q_aug'), 'Q_aug field must exist.');
            testCase.verifyTrue(isfield(augData, 'augInfo'), 'augInfo field must exist.');
        end
        
        function testAryuleBranch(testCase)
            % --- Test ARYULE branch ---
            testCase.general_config.augmentation_model_initializer.id = 'aryule';
            testCase.general_config.augmentation_model_initializer.model_params.model_order = 2;
            
            augData = get_augmentation_model(testCase.training_data, testCase.general_config, testCase.sampling_interval);
            
            testCase.verifyEqual(lower(augData.augmentation_type), 'aryule', 'Augmentation type should be aryule.');
            testCase.verifyTrue(isfield(augData, 'F_aug'), 'F_aug field must exist.');
            testCase.verifyTrue(isfield(augData, 'Q_aug'), 'Q_aug field must exist.');
            testCase.verifyTrue(isfield(augData, 'augInfo'), 'augInfo field must exist.');
        end
        
        function testKinematicBranch(testCase)
            % --- Test Kinematic branch ---
            testCase.general_config.augmentation_model_initializer.id = 'kinematic';
            testCase.general_config.augmentation_model_initializer.model_params.wiener_mdl_order = 3;
            testCase.general_config.augmentation_model_initializer.model_params.process_noise_variance = 0.1;
            
            augData = get_augmentation_model(testCase.training_data, testCase.general_config, testCase.sampling_interval);
            
            testCase.verifyEqual(lower(augData.augmentation_type), 'kinematic', 'Augmentation type should be kinematic.');
            testCase.verifyTrue(isfield(augData, 'F_aug'), 'F_aug field must exist.');
            testCase.verifyTrue(isfield(augData, 'Q_aug'), 'Q_aug field must exist.');
        end
        
        function testArimaBranch(testCase)
            % --- Test ARIMA branch ---
            testCase.general_config.augmentation_model_initializer.id = 'arima';
            testCase.general_config.augmentation_model_initializer.model_params.p = 1;
            testCase.general_config.augmentation_model_initializer.model_params.D = 1;
            testCase.general_config.augmentation_model_initializer.model_params.q = 1;
            
            % The local stubs for arima and estimate (defined below)
            % yield predictable ARIMA state-space matrices.
            augData = get_augmentation_model(testCase.training_data, testCase.general_config, testCase.sampling_interval);
            
            testCase.verifyEqual(lower(augData.augmentation_type), 'arima', 'Augmentation type should be arima.');
            testCase.verifyTrue(isfield(augData, 'F_aug'), 'F_aug field must exist.');
            testCase.verifyTrue(isfield(augData, 'Q_aug'), 'Q_aug field must exist.');
            testCase.verifyTrue(isfield(augData, 'H_aug'), 'H_aug field must exist for ARIMA.');
        end
        
        function testRbfBranch(testCase)
            % --- Test RBF branch (should error) ---
            testCase.general_config.augmentation_model_initializer.id = 'rbf';
            testCase.verifyError(@() get_augmentation_model(testCase.training_data, testCase.general_config, testCase.sampling_interval), ...
                'MATLAB:RBFUnavailable');
        end
        
        function testNoneBranch(testCase)
            % --- Test 'none' branch ---
            testCase.general_config.augmentation_model_initializer.id = 'none';
            augData = get_augmentation_model(testCase.training_data, testCase.general_config, testCase.sampling_interval);
            testCase.verifyEqual(lower(augData.augmentation_type), 'none', 'Augmentation type should be none.');
            % For 'none', we expect no extra augmentation matrices.
            testCase.verifyEqual(length(fieldnames(augData)), 1, 'Only augmentation_type field should be present in "none" branch.');
        end
        
        function testInvalidBranch(testCase)
            % --- Test an invalid augmentation model ---
            testCase.general_config.augmentation_model_initializer.id = 'invalid_model';
            testCase.verifyError(@() get_augmentation_model(testCase.training_data, testCase.general_config, testCase.sampling_interval), ...
                'MATLAB:UndefinedAugmentationModel');
        end
    end
end

%--------------------------------------------------------------------------
% Local stub functions for test_get_augmentation_model

function [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = arfit(training_data, model_order1, model_order2)
    % Dummy arfit stub
    intercept_vector = 0.1;
    var_coefficient_matrices = ones(2);
    var_covariance_matrices = eye(2);
end

function [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(~, ~)
    % Dummy construct_var_matrices stub
    F_var = eye(2)*2; 
    Q_var = eye(2)*3; 
    var_states_amount = 2; 
    var_model_order = 2;
end

function [ar_coefficients, ar_variance] = aryule(training_data, model_order)
    % Dummy aryule stub
    ar_coefficients = [1, 0.5, 0.5];
    ar_variance = 1;
end

function [F_wiener_aug, Q_wiener_aug] = get_discrete_wiener_model(varargin)
    % Dummy stub for kinematic branch.
    F_wiener_aug = eye(2)*5;
    Q_wiener_aug = eye(2)*6;
end

function mdl = arima(varargin)
    % Dummy stub for arima constructor. Return a dummy struct.
    mdl = struct();
end

function est_mdl = estimate(mdl, training_data)
    % Dummy estimate stub.
    % Simulate estimated AR, MA coefficients and variance.
    est_mdl.AR = {1}; %#ok<CCAT1>
    est_mdl.MA = {0.5}; %#ok<CCAT1>
    est_mdl.Variance = 1;
end

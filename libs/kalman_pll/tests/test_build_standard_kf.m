classdef test_build_standard_kf < matlab.unittest.TestCase
    properties
        F_los
        Q_los
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
            % Create a simple LOS dynamics example.
            testCase.F_los = eye(2);
            testCase.Q_los = eye(2);
            testCase.sampling_interval = 0.01;
            testCase.general_config.C_over_N0_array_dBHz = 35;
        end
    end

    methods(Test)
        function testNoAugmentation(testCase)
            % --- Test the 'none' branch ---
            aug_data.augmentation_type = 'none';
            
            [F, Q, H, R, W] = build_standard_kf(testCase.F_los, testCase.Q_los, aug_data, testCase.general_config, testCase.sampling_interval);
            
            % Expect LOS matrices to be used directly.
            testCase.verifyEqual(F, testCase.F_los, 'F should equal F_los when no augmentation is applied.');
            testCase.verifyEqual(Q, testCase.Q_los, 'Q should equal Q_los when no augmentation is applied.');
            
            % H is constructed as [1, zeros(1, n_los-1)]
            expectedH = [1, 0]; 
            testCase.verifyEqual(H, expectedH, 'H must be [1,0] for LOS only.');
            
            testCase.verifyEqual(W, zeros(2, 1), 'W must be a zero vector of length equal to n_los.');
        end
        
        function testArfitAugmentation(testCase)
            % --- Test ARFIT branch ---
            aug_data.augmentation_type = 'arfit';
            aug_data.F_aug = eye(3)*2;
            aug_data.Q_aug = eye(3)*3;
            aug_data.intercept = 0.5;
            aug_data.augInfo.states_amount = 1;
            aug_data.augInfo.model_order = 3;
            
            % When using ARFIT (or aryule) branch, build_standard_kf calls construct_kalman_matrices.
            % The stub (defined below) returns fixed matrices.
            [F, Q, H, R, W] = build_standard_kf(testCase.F_los, testCase.Q_los, aug_data, testCase.general_config, testCase.sampling_interval);
            
            % The stub returns:
            %   F = [10 0; 0 10]
            %   Q = [20 0; 0 20]
            %   H = [1 0]
            %   R = diag([40, 40])
            %   W = [50; 50]
            stubF = diag([1,1,2,2,2]);
            stubQ = diag([1,1,3,3,3]);
            stubH = [1, 0, 1, 0, 0];
            stubR = get_phase_variances(testCase.general_config.C_over_N0_array_dBHz, ...
                testCase.sampling_interval);
            stubW = [0, 0, 0.5, 0, 0].';
            
            testCase.verifyEqual(F, stubF, 'F must match the stubbed output from construct_kalman_matrices.');
            testCase.verifyEqual(Q, stubQ, 'Q must match the stubbed output.');
            testCase.verifyEqual(H, stubH, 'H must match the stubbed output.');
            testCase.verifyEqual(R, stubR, 'R must match the stubbed output.');
            testCase.verifyEqual(W, stubW, 'W must match the stubbed output.');
        end
        
        function testKinematicAugmentation(testCase)
            % --- Test kinematic branch ---
            aug_data.augmentation_type = 'kinematic';
            % Provide dummy augmentation matrices for the kinematic branch.
            aug_data.F_aug = eye(2)*5;
            aug_data.Q_aug = eye(2)*6;
            
            [F, Q, H, R, W] = build_standard_kf(testCase.F_los, testCase.Q_los, aug_data, testCase.general_config, testCase.sampling_interval);
            
            % Expected values: F and Q are block-diagonal combinations.
            expectedF = blkdiag(testCase.F_los, aug_data.F_aug);
            expectedQ = blkdiag(testCase.Q_los, aug_data.Q_aug);
            testCase.verifyEqual(F, expectedF, 'F must be block-diagonal of LOS and augmentation matrices.');
            testCase.verifyEqual(Q, expectedQ, 'Q must be block-diagonal of LOS and augmentation matrices.');
            
            n_los = size(testCase.F_los,1);
            n_aug = size(aug_data.F_aug,1);
            expectedH = [1, zeros(1, n_los-1), 1, zeros(1, n_aug-1)];
            testCase.verifyEqual(H, expectedH, 'H must be correctly constructed for kinematic augmentation.');
            testCase.verifyEqual(W, zeros(n_los+n_aug, 1), 'W must be zeros.');
        end
        
        function testArimaAugmentation(testCase)
            % --- Test ARIMA branch ---
            aug_data.augmentation_type = 'arima';
            aug_data.F_aug = eye(3)*7;
            aug_data.Q_aug = eye(3)*8;
            aug_data.H_aug = [0.5, 0.5, 0];
            
            [F, Q, H, R, W] = build_standard_kf(testCase.F_los, testCase.Q_los, aug_data, testCase.general_config, testCase.sampling_interval);
            
            expectedF = blkdiag(testCase.F_los, aug_data.F_aug);
            expectedQ = blkdiag(testCase.Q_los, aug_data.Q_aug);
            testCase.verifyEqual(F, expectedF, 'F must be block-diagonal for ARIMA augmentation.');
            testCase.verifyEqual(Q, expectedQ, 'Q must be block-diagonal for ARIMA augmentation.');
            
            expectedH = [1, zeros(1, size(testCase.F_los,1)-1), aug_data.H_aug];
            testCase.verifyEqual(H, expectedH, 'H must combine LOS H and ARIMA H.');
            testCase.verifyEqual(W, zeros(size(F,1),1), 'W must be a zero vector for ARIMA branch.');
        end
        
        function testInvalidAugmentationType(testCase)
            % --- Test for unsupported augmentation type ---
            aug_data.augmentation_type = 'unsupported';
            testCase.verifyError(@() build_standard_kf(testCase.F_los, testCase.Q_los, aug_data, testCase.general_config, testCase.sampling_interval), ...
                'MATLAB:UndefinedKFType');
        end
    end
end
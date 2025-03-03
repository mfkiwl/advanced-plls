classdef test_construct_kalman_matrices < matlab.unittest.TestCase
% test_construct_kalman_matrices
%
% This test suite verifies the behavior of the function construct_kalman_matrices,
% ensuring that it correctly constructs Kalman filter matrices while enforcing
% valid constraints on input parameters.
%
% Tests include:
%   - A valid input scenario.
%   - Verification of expected matrix dimensions and values.
%   - Error cases for invalid inputs, such as:
%       - Non-numeric, non-square, or non-positive semi-definite covariance matrices.
%       - Invalid intercept vector.
%       - Invalid process noise covariance matrices.

    methods(TestClassSetup)
        function addParentPath(testCase)
            % Add the parent directory to the path (if necessary).
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parent_dir);
            testCase.addTeardown(@() rmpath(parent_dir));
        end
    end

    methods(Test)
        function testValidInputs(testCase)
            % Test valid Kalman filter matrix construction.

            % Define valid LOS matrices (2x2)
            F_los = eye(2);
            Q_los = [0.1, 0.02; 0.02, 0.1];  % Symmetric and positive semi-definite

            % Define valid VAR matrices for VAR(2) with 2 states
            var_states_amount = 2;
            var_model_order = 2;
            F_var = [0.5, 0.1, -0.2, 0.3; 
                     0.2, 0.3,  0.1, 0.4];
            Q_var = blkdiag([0.05, 0.01; 0.01, 0.05], zeros(2)); % Symmetric and PSD

            % Valid intercept vector
            intercept_vector = [0.1; -0.1];

            % Other parameters
            C_over_N0_array_dBHz = 35;
            sampling_interval = 0.01;

            % Call the function
            [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, var_states_amount, var_model_order, C_over_N0_array_dBHz, sampling_interval);

            % Expected outputs
            expected_F = blkdiag(F_los, F_var);
            expected_Q = blkdiag(Q_los, Q_var);
            expected_H = [1, 0, 1, zeros(1, var_states_amount * var_model_order - 1)];
            expected_R = diag(compute_phase_variances(C_over_N0_array_dBHz, sampling_interval));
            expected_W = [zeros(size(F_los,1),1); intercept_vector; zeros(var_states_amount * (var_model_order - 1),1)];

            % Assertions
            testCase.verifyEqual(F, expected_F, 'AbsTol', 1e-10);
            testCase.verifyEqual(Q, expected_Q, 'AbsTol', 1e-10);
            testCase.verifyEqual(H, expected_H, 'AbsTol', 1e-10);
            testCase.verifyEqual(R, expected_R, 'AbsTol', 1e-10);
            testCase.verifyEqual(W, expected_W, 'AbsTol', 1e-10);
        end

        function testNonSymmetricQ_los(testCase)
            % Q_los must be symmetric; this should trigger an error.
            F_los = eye(2);
            Q_los = [0.1, 0.02; 0.03, 0.1]; % Non-symmetric
            F_var = eye(2);
            Q_var = eye(2);
            intercept_vector = [0.1; -0.1];

            testCase.verifyError(@() construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, 2, 2, 35, 0.01), ...
                'construct_kalman_matrices:QlosNotSymmetric');
        end

        function testNonSymmetricQ_var(testCase)
            % Q_var must be symmetric; this should trigger an error.
            F_los = eye(2);
            Q_los = eye(2);
            F_var = eye(2);
            Q_var = [0.05, 0.01; 0.02, 0.05]; % Non-symmetric
            intercept_vector = [0.1; -0.1];

            testCase.verifyError(@() construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, 2, 2, 35, 0.01), ...
                'construct_kalman_matrices:QvarNotSymmetric');
        end

        function testNegativeDefiniteQ_los(testCase)
            % Q_los must be positive semi-definite; this should trigger an error.
            F_los = eye(2);
            Q_los = [-1, 0; 0, -1]; % Negative definite
            F_var = eye(2);
            Q_var = eye(2);
            intercept_vector = [0.1; -0.1];

            testCase.verifyError(@() construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, 2, 2, 35, 0.01), ...
                'construct_kalman_matrices:QlosNotPositiveSemiDefinite');
        end

        function testNegativeDefiniteQ_var(testCase)
            % Q_var must be positive semi-definite; this should trigger an error.
            F_los = eye(2);
            Q_los = eye(2);
            F_var = eye(2);
            Q_var = [-1, 0; 0, -1]; % Negative definite
            intercept_vector = [0.1; -0.1];

            testCase.verifyError(@() construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, 2, 2, 35, 0.01), ...
                'construct_kalman_matrices:QvarNotPositiveSemiDefinite');
        end

        function testInvalidInterceptVector(testCase)
            % Intercept vector must be a column; this should trigger an error.
            F_los = eye(2);
            Q_los = eye(2);
            F_var = eye(2);
            Q_var = eye(2);
            intercept_vector = [0.1, -0.1]; % Row vector (invalid)

            testCase.verifyError(@() construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, 2, 2, 35, 0.01), ...
                'MATLAB:construct_kalman_matrices:expectedColumn');
        end
    end
end

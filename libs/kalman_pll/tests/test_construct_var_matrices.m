classdef test_construct_var_matrices < matlab.unittest.TestCase
% test_construct_var_matrices
%
% This test suite verifies the behavior of the function
% construct_var_matrices, which constructs augmented VAR model
% matrices from given coefficient and covariance matrices.
%
% The tests include:
%   - Valid inputs for a VAR(1) model.
%   - Valid inputs for a VAR(2) model.
%   - An invalid case where the number of columns in the coefficient
%     matrix is not an integer multiple of the number of states.
%   - Non-numeric coefficient matrix input.
%   - Non-square covariance matrix input.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    methods(TestClassSetup)
        function add_parent_path(test_case)
            % Add the parent directory path
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parent_dir);

            test_case.addTeardown(@() rmpath(parent_dir));
        end
    end

    methods(Test)
        function testVAR1(testCase)
            % Test a valid VAR(1) model.
            % For a VAR(1), the coefficient matrix is n x n.
            n = 2;
            p = 1;
            % Construct a 2x2 coefficient matrix.
            A = [0.5, 0.1; 0.2, 0.3];
            % Construct a 2x2 covariance matrix.
            sigma = [1, 0.2; 0.2, 0.5];
            
            [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(A, sigma);
            
            % For VAR(1), F_var should equal A, and Q_var equals sigma.
            testCase.verifyEqual(F_var, A, 'AbsTol', 1e-10);
            testCase.verifyEqual(Q_var, sigma, 'AbsTol', 1e-10);
            testCase.verifyEqual(var_states_amount, n);
            testCase.verifyEqual(var_model_order, p);
        end
        
        function testVAR2(testCase)
            % Test a valid VAR(2) model.
            % For a VAR(2) model with n = 2, the coefficient matrix should be 2 x (2*2) = 2 x 4.
            n = 2;
            p = 2;
            A = [0.5, 0.1, -0.2, 0.3;
                 0.2, 0.3,  0.1, 0.4];  % 2x4 matrix.
            sigma = [1, 0.2; 0.2, 0.5];  % 2x2 covariance matrix.
            
            [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(A, sigma);
            
            % Expected augmented F_var:
            % [ A ;
            %   [I_n, zeros(n)] ] i.e., [A; [eye(2), zeros(2,2)]]
            expected_F_var = [A; [eye(n), zeros(n,n)]];
            % Expected Q_var = blkdiag(sigma, zeros(n))
            expected_Q_var = blkdiag(sigma, zeros(n));
            
            testCase.verifyEqual(F_var, expected_F_var, 'AbsTol', 1e-10);
            testCase.verifyEqual(Q_var, expected_Q_var, 'AbsTol', 1e-10);
            testCase.verifyEqual(var_states_amount, n);
            testCase.verifyEqual(var_model_order, p);
        end
        
        function testInvalidDimensions(testCase)
            % Test with invalid dimensions for var_coefficient_matrices.
            % For n = 2, the number of columns must be a multiple of 2.
            % Here, we provide a 2x3 matrix, which is invalid.
            A = [0.5, 0.1, -0.2; 0.3, 0.4, 0.5];  % 2x3 matrix (3 is not a multiple of 2)
            cov = [1, 0.2; 0.2, 0.5];  % 2x2 matrix.
            
            testCase.verifyError(@() construct_var_matrices(A, cov), ...
                'construct_var_matrices:InvalidDimensions', ...
                'Expected an error when the number of columns in var_coefficient_matrices is not a multiple of the number of rows.');
        end

        function testNonSquareCovariance(testCase)
            % Test with a non-square covariance matrix.
            % For a valid covariance matrix with n = 2, it must be 2x2.
            % Here, we provide a 2x3 matrix.
            A = [0.5, 0.1; 0.2, 0.3];  % 2x2 (VAR(1))
            cov = [1, 0.2, 0.3; 0.2, 0.5, 0.1];  % 2x3 matrix, invalid.
            
            testCase.verifyError(@() construct_var_matrices(A, cov), ...
                'MATLAB:construct_var_matrices:expectedSquare');
        end

        function testNonNumericCoefficient(testCase)
            % Test with a non-numeric coefficient matrix.
            A = 'not a numeric matrix';
            sigma = [1, 0.2; 0.2, 0.5];
            
            % Expect an error from validateattributes.
            testCase.verifyError(@() construct_var_matrices(A, sigma), ...
                'MATLAB:construct_var_matrices:invalidType');
        end
    end
end

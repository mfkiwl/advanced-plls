classdef test_get_discrete_wiener_model < matlab.unittest.TestCase
    % test_get_discrete_wiener_model
    %
    % Syntax:
    %   results = runtests('test_get_discrete_wiener_model')
    %
    % Description:
    %   Unit tests for get_discrete_wiener_model, which generates the 
    %   state transition and noise covariance matrices for LOS dynamics.
    %   Validates correct handling of inputs, output sizes, and error 
    %   handling on invalid data.
    %
    % Example:
    %   % Run the test suite:
    %   results = runtests('test_get_discrete_wiener_model');
    %   disp(results);
    %
    % Author:
    %   Rodrigo de Lima Florindo
    %   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com

    methods(TestClassSetup)
        function add_parent_path(test_case)
            % add_parent_path - Add the parent directory containing get_discrete_wiener_model.m
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parent_dir);
            test_case.addTeardown(@() rmpath(parent_dir));
        end
    end

    methods(Test)
        function test_valid_input_simple_case(test_case)
            % Test with a small valid scenario: L=2, M=2, etc.
            L = 2; 
            M = 2;
            sampling_interval = 0.01;
            sigma_array = [1e-3, 1e-4, 5e-4];
            delta_array = [1, 0.9];

            [stm, covm] = get_discrete_wiener_model(L, M, sampling_interval, sigma_array, delta_array);

            % Check sizes
            expected_dim = L + M - 1; % 2+2-1=3
            test_case.verifySize(stm, [expected_dim, expected_dim], ...
                'State transition matrix has incorrect size.');
            test_case.verifySize(covm, [expected_dim, expected_dim], ...
                'Covariance matrix has incorrect size.');
        end

        function test_edge_case_l1(test_case)
            % If L=1, the function might produce a degenerate scenario. 
            % Check if it handles that or triggers an error gracefully.
            L = 1;
            M = 2;
            sampling_interval = 0.01;
            sigma_array = [1e-3, 5e-4]; % L+M-1=1+2-1=2
            delta_array = 1;           % L=1

            [stm, covm] = get_discrete_wiener_model(L, M, sampling_interval, sigma_array, delta_array);

            expected_dim = L + M - 1; % 1+2-1=2
            test_case.verifySize(stm, [expected_dim, expected_dim], ...
                'STM dimension mismatch for L=1.');
            test_case.verifySize(covm, [expected_dim, expected_dim], ...
                'Covariance dimension mismatch for L=1.');
        end
        
        function test_invalid_l(test_case)
            % Non-integer or negative L should throw an error
            invalid_l = -1;
            M = 2;
            sampling_interval = 0.01;
            sigma_array = [1e-3, 5e-4];
            delta_array = 1;
            
            test_case.verifyError(@() get_discrete_wiener_model(invalid_l, M, sampling_interval, sigma_array, delta_array), ...
                'MATLAB:get_discrete_wiener_model:expectedPositive', ...
                'Expected an error for negative L.');
        end

        function test_invalid_sigma_length(test_case)
            % If sigma_array length != L+M-1, must throw an error
            L = 2; M = 2;
            sampling_interval = 0.01;
            sigma_array = 1e-3; % Should have 3 elements, not 1
            delta_array = [1, 0.9];
            
            test_case.verifyError(@() get_discrete_wiener_model(L, M, sampling_interval, sigma_array, delta_array), ...
                'MATLAB:get_discrete_wiener_model:incorrectNumel', ...
                'Expected an error due to mismatch in sigma_array length.');
        end
        
        function test_invalid_delta_size(test_case)
            % If delta_array length != L, must throw an error
            L = 2; M = 2;
            sampling_interval = 0.01;
            sigma_array = [1e-3, 1e-4, 5e-4]; % correct length: L+M-1=3
            delta_array = [1, 0.9, 0.85]; % Should have 2 elements, not 3

            test_case.verifyError(@() get_discrete_wiener_model(L, M, sampling_interval, sigma_array, delta_array), ...
                'MATLAB:get_discrete_wiener_model:incorrectNumel', ...
                'Expected an error for invalid delta_array size.');
        end
    end
end

classdef test_struct_to_nv_pairs < matlab.unittest.TestCase
% test_struct_to_nv_pairs
% Unit tests for the struct_to_nv_pairs function.
%
% Syntax:
%   results = runtests('test_struct_to_nv_pairs')
%
% Description:
%   This test suite verifies that a nonempty struct is correctly converted
%   into an alternating cell array of name-value pairs. It also tests that
%   invalid inputs (empty structs or non-structs) trigger appropriate errors.
%
% Example:
%   % Run the test suite:
%   results = runtests('test_struct_to_nv_pairs');
%   disp(results);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    methods(TestClassSetup)
        function add_parent_path(test_case)
            % add_parent_path - Add the parent directory containing struct_to_nv_pairs.
            % Since this test file is in "kalman_pll/tests", the parent directory is one level up.
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parent_dir);
            test_case.addTeardown(@() rmpath(parent_dir));
        end
    end

    methods(Test)
        function test_valid_conversion(test_case)
            % Verify that a simple struct is correctly converted to name-value pairs.
            input_struct = struct('a', 1, 'b', 2, 'c', 3);
            nv_pairs = struct_to_nv_pairs(input_struct);
            expected = {'a', 1, 'b', 2, 'c', 3};
            test_case.verifyEqual(nv_pairs, expected, 'The name-value pair output is incorrect.');
        end
        
        function test_order_and_length(test_case)
            % Verify that the output length equals twice the number of fields and order is preserved.
            input_struct = struct('x', 10, 'y', 20);
            nv_pairs = struct_to_nv_pairs(input_struct);
            field_names = fieldnames(input_struct);
            test_case.verifyEqual(numel(nv_pairs), 2 * numel(field_names), ...
                'The length of the name-value array is not correct.');
            for k = 1:numel(field_names)
                test_case.verifyEqual(nv_pairs{2*k-1}, field_names{k}, ...
                    sprintf('Field name at position %d is incorrect.', 2*k-1));
                test_case.verifyEqual(nv_pairs{2*k}, input_struct.(field_names{k}), ...
                    sprintf('Field value for %s is incorrect.', field_names{k}));
            end
        end
        
        function test_empty_struct_error(test_case)
            % Verify that an empty struct (with no fields) triggers an error.
            empty_struct = struct();  % This returns a 1x0 struct.
            test_case.verifyError(@() struct_to_nv_pairs(empty_struct), 'struct_to_nv_pairs:EmptyStruct');
        end
        
        function test_nonstruct_input_error(test_case)
            % Verify that a non-struct input triggers an error.
            non_struct = 42;
            % Update the expected error ID to match the one thrown by validateattributes.
            test_case.verifyError(@() struct_to_nv_pairs(non_struct), 'MATLAB:struct_to_nv_pairs:invalidType');
        end
    end
end

function name_value = struct_to_nv_pairs(s)
% struct_to_nv_pairs
%
% Syntax:
%   name_value = struct_to_nv_pairs(s)
%
% Description:
%   Converts the input struct into an alternating cell array of name-value pairs.
%   The output is organized as: {name1, value1, name2, value2, ...}.
%
% Inputs:
%   s - A nonempty struct with at least one field.
%
% Outputs:
%   name_value - A 1xN cell array containing the field names and corresponding values.
%
% Notes:
%   - The function validates that the input is a nonempty struct.
%   - If the struct has no fields, an error is thrown.
%
% Example:
%   % Convert a struct to name-value pairs:
%   s = struct('a', 1, 'b', 2, 'c', 3);
%   nv = struct_to_nv_pairs(s);
%   % nv will be: {'a', 1, 'b', 2, 'c', 3}
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Validate that s is a nonempty struct.
    validateattributes(s, {'struct'}, {'nonempty'}, mfilename, 's');
    % Check if the struct has at least one field.
    if isempty(fieldnames(s))
       error('struct_to_nv_pairs:EmptyStruct', 'Input struct must have at least one field.');
    end
    names = fieldnames(s);
    values = struct2cell(s);
    name_value = reshape([names, values]', 1, []);
end

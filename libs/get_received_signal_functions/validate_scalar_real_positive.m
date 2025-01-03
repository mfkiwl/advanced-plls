% validate_scalar_real_positive
% Validates whether an input is a positive real numeric scalar.
%
% Syntax:
%   validate_scalar_real_positive(value, function_name, name)
%
% Description:
%   This utility function checks if the given input `value` is a positive real
%   numeric scalar. If the validation fails, it raises an error with a specific
%   message and error identifier using the provided `function_name` and `name`.
%
% Inputs:
%   - value: The value to be validated.
%   - function_name: A string representing the name of the calling function, used
%                    in the error identifier.
%   - name: A string describing the name of the variable being validated, used
%           in the error message.
%
% Errors:
%   - Throws an error with an identifier of the form `<function_name>:InvalidInput`
%     if `value` is not numeric, not a scalar, or less than or equal to zero.
%
% Example:
%   validate_scalar_real_positive(3.14, 'example_function', 'example_variable');
%   % Passes since 3.14 is a positive scalar.
%
%   validate_scalar_real_positive(-5, 'example_function', 'example_variable');
%   % Raises an error: example_function:InvalidInput
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

function validate_scalar_real_positive(value, function_name, name)
    if ~isnumeric(value) || ~isscalar(value) || value <= 0
        error(sprintf('%s:InvalidInput', function_name), ...
              'The input "%s" must be a positive scalar. Received: %g.', name, num2str(value));
    end
end
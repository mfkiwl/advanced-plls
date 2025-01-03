function validate_scalar_numeric(value, function_name, name)
    if ~isnumeric(value) || ~isscalar(value)
        error(sprintf('%s:InvalidInput', function_name), ...
              'The input "%s" must be a numeric scalar.', name);
    end
end
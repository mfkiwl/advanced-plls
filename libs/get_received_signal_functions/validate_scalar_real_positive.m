function validate_scalar_real_positive(value, function_name, name)
    if ~isnumeric(value) || ~isscalar(value) || value <= 0
        error(sprintf('%s:InvalidInput', function_name), ...
              'The input "%s" must be a positive scalar. Received: %g.', name, num2str(value));
    end
end
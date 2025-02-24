function validate_config(config)
% validate_config
% Validates the input `config` struct for the Kalman PLL configuration.
%
% Throws an error if required fields are missing or if any field does not meet
% the expected data type and value constraints.

    % Required fields in the top-level config struct
    required_fields = {'discrete_wiener_model_config', ...
                       'scintillation_training_data_config', ...
                       'var_minimum_order', ...
                       'var_maximum_order', ...
                       'C_over_N0_array_dBHz', ...
                       'initial_states_distributions_boundaries', ...
                       'real_doppler_profile', ...
                       'is_use_cached_settings', ...
                       'is_generate_random_initial_estimates'};

    % Check for missing fields in the top-level config.
    missing_fields = setdiff(required_fields, fieldnames(config));
    if ~isempty(missing_fields)
        error('validate_config:MissingFields', ...
              'The following required fields are missing: %s', strjoin(missing_fields, ', '));
    end

    % Validate attributes for each top-level field.
    try
        validateattributes(config.discrete_wiener_model_config, {'cell'}, {'nonempty'}, 'validate_config', 'discrete_wiener_model_config');
        % Now scintillation_training_data_config is expected to be a nonempty struct.
        validateattributes(config.scintillation_training_data_config, {'struct'}, {'nonempty'}, 'validate_config', 'scintillation_training_data_config');
        validateattributes(config.var_minimum_order, {'numeric'}, {'scalar', 'integer', 'positive', 'nonnan'}, 'validate_config', 'var_minimum_order');
        validateattributes(config.var_maximum_order, {'numeric'}, {'scalar', 'integer', '>=', config.var_minimum_order}, 'validate_config', 'var_maximum_order');
        validateattributes(config.C_over_N0_array_dBHz, {'numeric'}, {'vector', 'real', 'finite', 'nonnan', 'nonempty'}, 'validate_config', 'C_over_N0_array_dBHz');
        validateattributes(config.initial_states_distributions_boundaries, {'cell'}, {'nonempty'}, 'validate_config', 'initial_states_distributions_boundaries');
        validateattributes(config.real_doppler_profile, {'numeric'}, {'vector', 'real', 'finite', 'nonnan'}, 'validate_config', 'real_doppler_profile');
        validateattributes(config.is_use_cached_settings, {'logical'}, {'scalar'}, 'validate_config', 'is_use_cached_settings');
        validateattributes(config.is_generate_random_initial_estimates, {'logical'}, {'scalar'}, 'validate_config', 'is_generate_random_initial_estimates');
    catch ME
        rethrow(ME);
    end

    % Additional consistency checks on top-level config.
    if config.var_minimum_order > config.var_maximum_order
        error('validate_config:InvalidOrderRange', ...
              '`var_minimum_order` (%d) cannot be greater than `var_maximum_order` (%d).', ...
              config.var_minimum_order, config.var_maximum_order);
    end

    % Validate scintillation_training_data_config.
    tpStruct = config.scintillation_training_data_config;
    if ~isfield(tpStruct, 'scintillation_model')
        error('validate_config:InvalidScintDataFormat', ...
              'scintillation_training_data_config must contain a ''scintillation_model'' field.');
    end
    model = tpStruct.scintillation_model;
    validModels = {'CSM', 'TPPSM', 'none'};
    if ~ismember(upper(model), upper(validModels))
        error('validate_config:InvalidScintDataFormat', ...
              'scintillation_model must be one of: %s', strjoin(validModels, ', '));
    end

    switch upper(model)
        case 'TPPSM'
            % For TPPSM, required fields: scenario, simulation_time, sampling_interval,
            % and is_refractive_effects_removed.
            requiredTPPSMFields = {'scenario', 'simulation_time', 'sampling_interval', 'is_refractive_effects_removed'};
            missingTPPSMFields = setdiff(requiredTPPSMFields, fieldnames(tpStruct));
            if ~isempty(missingTPPSMFields)
                error('validate_config:InvalidScintDataFormat', ...
                      'For TPPSM, scintillation_training_data_config must contain the fields: %s', strjoin(requiredTPPSMFields, ', '));
            end
            % Validate scenario.
            validateattributes(tpStruct.scenario, {'char', 'string'}, {'nonempty'}, 'validate_config', 'scintillation_training_data_config scenario');
            validScenarios = {'Weak', 'Moderate', 'Severe'};
            if ~ismember(char(tpStruct.scenario), validScenarios)
                error('validate_config:InvalidScintDataFormat', ...
                      'For TPPSM, scenario must be one of: %s', strjoin(validScenarios, ', '));
            end
            % Validate simulation_time and sampling_interval.
            validateattributes(tpStruct.simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, ...
                'validate_config', 'scintillation_training_data_config simulation_time (TPPSM)');
            validateattributes(tpStruct.sampling_interval, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, ...
                'validate_config', 'scintillation_training_data_config sampling_interval (TPPSM)');
            % Validate is_refractive_effects_removed.
            validateattributes(tpStruct.is_refractive_effects_removed, {'logical'}, {'scalar'}, ...
                'validate_config', 'scintillation_training_data_config is_refractive_effects_removed (TPPSM)');
            
        case 'CSM'
            % For CSM, required fields: S4, tau0, simulation_time, sampling_interval.
            requiredCSMFields = {'S4', 'tau0', 'simulation_time', 'sampling_interval'};
            missingCSMFields = setdiff(requiredCSMFields, fieldnames(tpStruct));
            if ~isempty(missingCSMFields)
                error('validate_config:InvalidScintDataFormat', ...
                      'For CSM, scintillation_training_data_config must contain the fields: %s', strjoin(requiredCSMFields, ', '));
            end
            validateattributes(tpStruct.S4, {'numeric'}, {'scalar', 'real', 'finite', 'nonnan', '>=', 0, '<=', 1}, ...
                'validate_config', 'scintillation_training_data_config S4');
            validateattributes(tpStruct.tau0, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, ...
                'validate_config', 'scintillation_training_data_config tau0');
            validateattributes(tpStruct.simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, ...
                'validate_config', 'scintillation_training_data_config simulation_time (CSM)');
            validateattributes(tpStruct.sampling_interval, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, ...
                'validate_config', 'scintillation_training_data_config sampling_interval (CSM)');
            
        case 'NONE'
            % For 'none', you might not need additional fields.
            % (Optionally, you can enforce that no extra fields are present.)
            % Here we assume no extra validation is needed.
    end

    % Check that initial_states_distributions_boundaries and real_doppler_profile
    % have the same number of elements.
    if numel(config.initial_states_distributions_boundaries) ~= numel(config.real_doppler_profile)
        error('validate_config:BoundaryProfileMismatch', ...
              '`initial_states_distributions_boundaries` must have the same number of elements as `real_doppler_profile`. Found %d boundaries and %d profile elements.', ...
              numel(config.initial_states_distributions_boundaries), numel(config.real_doppler_profile));
    end

    % Check that each boundary is a valid range.
    for i = 1:numel(config.initial_states_distributions_boundaries)
        boundary = config.initial_states_distributions_boundaries{i};
        if ~isnumeric(boundary) || numel(boundary) ~= 2 || boundary(1) >= boundary(2)
            error('validate_config:InvalidBoundaryRange', ...
                  'Each element of `initial_states_distributions_boundaries` must be a numeric vector containing a valid range [a, b] where a < b.');
        end
    end
end

function kalman_pll_config = get_kalman_pll_config(config)
    % get_kalman_pll_config
    % Generates or retrieves Kalman filter parameters based on the provided configuration.
    %
    % Syntax:
    %   kalman_pll_config = get_kalman_pll_config(config)
    %
    % Description:
    %   This function computes or retrieves the Kalman filter parameters, including 
    %   state-space matrices and noise covariance matrices. It supports caching to 
    %   reuse previously computed parameters for efficient performance. If caching is 
    %   not enabled or the cache is invalid, it computes new parameters and updates the cache.
    %
    % Inputs:
    %   config - Struct containing all configuration details with the following fields:
    %       - discrete_wiener_model_config: Cell array for LOS dynamics parameters.
    %         Example: {1,3,0.01,[0,0,1],1}, where:
    %           * 1  - Number of frequency bands (L = 1),
    %           * 3  - Third-order LOS dynamics (M = 3),
    %           * 0.01 - Sampling time (in seconds),
    %           * [0,0,1] - Sigma array for noise levels,
    %           * 1  - Ratio of its frequency to a reference frequency (delta_array).
    %
    %       - scintillation_training_data_config: Cell array for scintillation model parameters.
    %         Example: {0.8, 0.7, 600, 0.01}, where:
    %           * 0.8  - S4 index,
    %           * 0.7  - Tau0 (decorrelation time, in seconds),
    %           * 600  - Total simulation time (in seconds),
    %           * 0.01 - Sampling time (in seconds).
    %
    %       - var_minimum_order: Minimum VAR (Vector Autoregressive) model order.
    %       - var_maximum_order: Maximum VAR model order.
    %       - C_over_N0_array_dBHz: Array representing the average C/N0 values for each 
    %         frequency band (in dB-Hz).
    %       - scint_model: Selected scintillation model ('CSM', 'MFPSM', or 'none').
    %       - is_refractive_effects_removed: Boolean flag indicating whether refractive 
    %         effects are removed for MFPSM.
    %       - is_use_cached_model: Boolean flag indicating whether cached configurations 
    %         should be used.
    %
    % Outputs:
    %   kalman_pll_config - Struct containing the computed Kalman filter parameters, with the following fields:
    %       * F  - Full state transition matrix.
    %       * Q  - Full process noise covariance matrix.
    %       * H  - Measurement matrix.
    %       * R  - Measurement noise covariance matrix.
    %       * F_los, Q_los - LOS dynamics matrices.
    %       * F_var, Q_var - VAR model matrices.
    %       * intercept_vector - VAR model intercept vector.
    %       * var_states_amount - Number of VAR model states.
    %       * var_model_order - Order of the VAR model.
    %
    % Notes:
    %   - Caching is used to improve performance by avoiding repeated computations.
    %   - Sampling intervals for LOS dynamics and scintillation models must match.
    %
    % Examples:
    %   % Example configuration:
    %   config = struct( ...
    %       'discrete_wiener_model_config', {1,3,0.01,[0,0,1],1}, ...
    %       'scintillation_training_data_config', {0.8, 0.7, 300, 0.01}, ...
    %       'var_minimum_order', 1, ...
    %       'var_maximum_order', 6, ...
    %       'C_over_N0_array_dBHz', 35, ...
    %       'scint_model', 'CSM', ...
    %       'is_refractive_effects_removed', false, ...
    %       'is_use_cached_model', false);
    %   kalman_pll_config = get_kalman_pll_config(config);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Cache file path
    cache_file = 'kalman_pll_cache.mat';

    % Validate required fields in config
    validate_config(config);

    % Extract sampling_interval from configurations for validation
    sampling_interval_los = config.discrete_wiener_model_config{1,3};
    sampling_interval_scint = config.scintillation_training_data_config{1,4};
    sampling_interval = validate_sampling_interval(sampling_interval_los, sampling_interval_scint);

    % Load or initialize kalman_pll_config
    [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file);

    % Handle computation or retrieval of Kalman parameters
    if is_cache_used
        fprintf('Using cached Kalman parameters for %s.\n', config.scint_model);
    else
        fprintf('Computing Kalman parameters for %s.\n', config.scint_model);

        % Compute model-specific parameters
        [F_los, Q_los] = get_discrete_wiener_model(config.discrete_wiener_model_config{:});

        if strcmp(config.scint_model, 'none')
            [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = deal([], [], F_los, Q_los, ...
                [1, zeros(1, size(F_los, 1) - 1)], diag(compute_phase_variances(config.C_over_N0_array_dBHz, sampling_interval)), ...
                [], [], []);
        else
            [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = compute_model_parameters( ...
                config.scint_model, config.scintillation_training_data_config, ...
                config.var_minimum_order, config.var_maximum_order, ...
                config.C_over_N0_array_dBHz, sampling_interval, F_los, ...
                Q_los, config.is_refractive_effects_removed);
        end

        % Store results in kalman_pll_config struct
        kalman_pll_config.(config.scint_model) = struct('F_los', F_los, 'Q_los', Q_los, ...
            'F_var', F_var, 'Q_var', Q_var, 'F', F, 'Q', Q, ...
            'H', H, 'R', R, 'intercept_vector', intercept_vector, 'var_model_order', ...
            var_model_order, 'var_states_amount', var_states_amount);

        % Save updated kalman_pll_config to cache
        save(cache_file, 'kalman_pll_config');
    end
end

function validate_config(config)
    % Ensure required fields are present in the configuration struct
    required_fields = {'discrete_wiener_model_config', ...
                       'scintillation_training_data_config', ...
                       'var_minimum_order', 'var_maximum_order', ...
                       'C_over_N0_array_dBHz', 'scint_model', ...
                       'is_refractive_effects_removed', ...
                       'is_use_cached_model'};
    missing_fields = setdiff(required_fields, fieldnames(config));
    if ~isempty(missing_fields)
        error('get_kalman_pll_config:MissingFields', ...
              'The following required fields are missing: %s', strjoin(missing_fields, ', '));
    end
end

function sampling_interval = validate_sampling_interval(sampling_interval_los, sampling_interval_scint)
    % validate_sampling_interval
    % Validates consistency of sampling intervals between LOS dynamics and scintillation models.
    %
    % Syntax:
    %   sampling_interval = validate_sampling_interval(sampling_interval_los, sampling_interval_scint)
    %
    % Description:
    %   This function ensures that the sampling intervals for LOS dynamics and 
    %   scintillation models match. If they differ, an error is raised.
    %
    % Inputs:
    %   sampling_interval_los   - Sampling interval for LOS dynamics (seconds).
    %   sampling_interval_scint - Sampling interval for scintillation model (seconds).
    %
    % Outputs:
    %   sampling_interval - The validated and consistent sampling interval.
    %
    % Notes:
    %   - The function assumes both intervals are scalars and checks for equality.
    %
    % Examples:
    %   sampling_interval = validate_sampling_interval(0.01, 0.01);
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Validate consistency of sampling_interval between los and scintillation kalman_pll_config
    if sampling_interval_los ~= sampling_interval_scint
        error(['Inconsistent sampling_interval values: los (%f) vs. ' ...
               'scintillation model (%f).'], sampling_interval_los, ...
               sampling_interval_scint);
    end
    sampling_interval = sampling_interval_los; % Return consistent value
end

function [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file)
    % load_or_initialize_models
    % Loads cached Kalman filter configuration or initializes a new one.
    %
    % Syntax:
    %   [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file)
    %
    % Description:
    %   This function attempts to load cached Kalman filter configurations 
    %   from a specified file. If the cache file does not exist or is invalid, 
    %   it initializes a new configuration.
    %
    % Inputs:
    %   config     - Struct containing user configuration.
    %   cache_file - Path to the cache file for storing Kalman filter configurations.
    %
    % Outputs:
    %   kalman_pll_config - The Kalman filter configuration struct.
    %   is_cache_used     - Boolean flag indicating whether cached values were used.
    %
    % Notes:
    %   - If the cache file exists but does not contain the required configuration 
    %     for the specified scintillation model, a new configuration is computed.
    %
    % Examples:
    %   [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, 'kalman_pll_cache.mat');
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Initialize the output flag
    is_cache_used = false;

    % Check if the cache file exists
    if isfile(cache_file)
        % Load existing cached kalman_pll_config
        fprintf('Cache file found. Loading cached kalman_pll_config.\n');
        load(cache_file, 'kalman_pll_config');

        % Check if the specific scintillation model in the cache is non-empty
        if isfield(kalman_pll_config, config.scint_model) && ...
           ~isempty(fieldnames(kalman_pll_config.(config.scint_model)))
            if config.is_use_cached_model
                fprintf('Using cached values for %s.\n', config.scint_model);
                is_cache_used = true; % Indicate cache was used
                return; % Use cached values
            else
                fprintf('Recomputing values for %s and updating cache.\n', config.scint_model);
            end
        else
            fprintf('Cache found but missing %s. Computing and caching new values.\n', config.scint_model);
        end
    else
        fprintf('No cache file found. Initializing kalman_pll_config.\n');
        % Initialize the kalman_pll_config struct
        kalman_pll_config = struct('CSM', struct(), 'MFPSM', struct(), 'none', struct());
    end

    % Always save to cache after initialization or recomputation
    if ~is_cache_used
        save(cache_file, 'kalman_pll_config');
        fprintf('Updated cache saved for %s.\n', config.scint_model);
    end
end

function [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = compute_model_parameters( ...
    scint_model, scintillation_training_data_config, var_minimum_order, ...
    var_maximum_order, C_over_N0_array_dBHz, sampling_interval, ...
    F_los, Q_los, is_refractive_effects_removed)
    % compute_model_parameters
    % Computes Kalman filter matrices for a specified scintillation model.
    %
    % Syntax:
    %   [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = ...
    %       compute_model_parameters(scint_model, scintillation_training_data_config, ...
    %       var_minimum_order, var_maximum_order, C_over_N0_array_dBHz, ...
    %       sampling_interval, F_los, Q_los, is_refractive_effects_removed)
    %
    % Description:
    %   This function computes the state-space matrices (F, Q, H, R) for the 
    %   Kalman filter based on the selected scintillation model (CSM or MFPSM).
    %   It also derives VAR model parameters.
    %
    % Inputs:
    %   scint_model                      - Selected scintillation model ('CSM', 'MFPSM').
    %   scintillation_training_data_config - Cell array for scintillation model parameters.
    %   var_minimum_order                - Minimum VAR model order.
    %   var_maximum_order                - Maximum VAR model order.
    %   C_over_N0_array_dBHz             - Array of C/N0 values in dB-Hz.
    %   sampling_interval                - Sampling interval for computations.
    %   F_los, Q_los                     - LOS dynamics matrices.
    %   is_refractive_effects_removed    - Boolean flag for MFPSM refractive effects.
    %
    % Outputs:
    %   F_var, Q_var     - State transition and covariance matrices for the VAR model.
    %   F, Q, H, R       - Full Kalman filter matrices.
    %   intercept_vector - Intercept vector for the VAR model.
    %   var_states_amount- Number of states in the VAR model.
    %   var_model_order  - Order of the VAR model.
    %
    % Notes:
    %   - The function handles training data preprocessing and VAR model fitting.
    %
    % Examples:
    %   [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = ...
    %       compute_model_parameters('CSM', scintillation_training_data_config, 1, 6, 35, 0.01, F_los, Q_los, false);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    switch scint_model
        case 'CSM'
            scint_complex_field = get_csm_data(scintillation_training_data_config{:});
            training_data = angle(scint_complex_field);

        case 'MFPSM'
            [scint_complex_field, ps_realization] = ...
                get_mfpsm_data(scintillation_training_data_config{:});
            if is_refractive_effects_removed
                warning('Removing refractive effects from MFPSM data.');
                scint_complex_field = scint_complex_field .* ...
                                      exp(-1j * ps_realization);
            end
            training_data = angle(scint_complex_field);

        otherwise
            error('Unsupported scintillation model: %s', scint_model);
    end

    % Fit var model
    [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = ...
        arfit(training_data, var_minimum_order, var_maximum_order);

    % Construct F_var and Q_var
    var_states_amount = size(var_coefficient_matrices, 1);
    var_model_order = size(var_coefficient_matrices, 2) / var_states_amount;

    F_var = [var_coefficient_matrices; [eye(var_states_amount * ...
        (var_model_order - 1)), zeros(var_states_amount * ...
        (var_model_order - 1), var_states_amount)]];
    Q_var = blkdiag(var_covariance_matrices, zeros(var_states_amount * ...
        (var_model_order - 1)));

    % Construct full Kalman matrices
    F = blkdiag(F_los, F_var);
    Q = blkdiag(Q_los, Q_var);
    H = [1, zeros(1, size(F_los, 1) - 1), 1, ...
         zeros(1, var_states_amount * var_model_order - 1)];
    R = diag(compute_phase_variances(C_over_N0_array_dBHz, sampling_interval));
end

function sigma2_array = compute_phase_variances(C_over_N0_array_dBHz, sampling_interval)
    % Compute phase variances based on C/N0
    C_over_N0_array_linear = 10.^(C_over_N0_array_dBHz ./ 10);
    sigma2_array = (1 ./ (2 * C_over_N0_array_linear * sampling_interval)) ...
        .* (1 + 1 ./ (2 * C_over_N0_array_linear * sampling_interval));
end

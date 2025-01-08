function [kalman_pll_config, initial_estimates] = get_kalman_pll_config(config)
    % get_kalman_pll_config
    % Generates or retrieves Kalman filter settings and initial estimates,
    % based on the provided configuration parameters.
    %
    % Syntax:
    %   kalman_pll_config = get_kalman_pll_config(config)
    %
    % Description:
    %   This function computes or retrieves the Kalman filter settings, including 
    %   state-space matrices and noise covariance matrices and its initial estimates. 
    %   It supports caching to reuse previously computed parameters for efficient 
    %   performance. If caching is not enabled or the cache is invalid, 
    %   it computes new parameters and updates the cache.
    %
    % Inputs:
    %   config - Struct containing all configuration details with the following fields:
    %       - discrete_wiener_model_config: Cell array for LOS dynamics parameters to 
    %                                       be used by `get_los_phase` function.
    %         Example: {1,3,0.01,[0,0,1],1}, where:
    %           * 1  - Number of frequency bands (`L`),
    %           * 3  - Third-order LOS dynamics (`M`),
    %           * 0.01 - Sampling time (`sampling_time`),
    %           * [0,0,1] - Sigma array for noise levels (`sigma_array`),
    %           * 1  - Ratio of its frequency to a reference frequency (`delta_array`).
    %
    %       - scintillation_training_data_config: Cell array for scintillation model parameters.
    %         Example: {0.8, 0.7, 600, 0.01}, where:
    %           * 0.8  - S4 index (`S4`),
    %           * 0.7  - τ₀ (`tau0`),
    %           * 600  - Total simulation time (`simulation_time`),
    %           * 0.01 - Sampling time (`sampling_time`).
    %
    %       - var_minimum_order: Minimum VAR (Vector Autoregressive) model
    %          order. The `arfit` function automatically estimates the vector
    %          autoregressive model order that lies withing a minimum and 
    %          maximum orders provided by the user that minimizes the 
    %          Schwarz Bayesian Criterion (SBC). It is important to comment
    %          that for time series with extremely large samples amount, the
    %          model order estimated by this function is generally close to
    %          the maximum one, since the SBC only starts to increase again at
    %          larger orders in this case.
    %       - var_maximum_order: Maximum VAR model order.
    %       - C_over_N0_array_dBHz: Array representing the average C/N0 values for each 
    %          frequency band (in dB-Hz).
    %       - training_scint_model: Selected scintillation model ('CSM', 'MFPSM', or 'none').
    %       - initial_states_distributions_boundaries: Uniform
    %          distributions boundaries for generating the initial state
    %          estimates for the Doppler profile used by the Kalman Filter.
    %       - real_doppler_profile: Real Doppler profile used to simulate the
    %          synthetic line-of-sight dynamics.
    %       - is_refractive_effects_removed: Boolean flag indicating whether refractive 
    %          effects are removed for MFPSM.
    %       - is_use_cached_settings: Boolean flag indicating whether cached configurations 
    %          should be used.
    %       - is_generate_random_initial_estimates: Boolean flag indicating
    %          if the inital estimates will be perfect or randomly generated
    %          according to the parsed
    %          `initial_states_distributions_boundaries` array.
    %
    % Outputs:
    %   kalman_pll_config - Struct containing the computed Kalman filter 
    %                       settings, with the following fields:
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
    %       'training_scint_model', 'CSM', ...
    %       'initial_states_distributions_boundaries', {[-pi,pi],[-5,5],[-0.1,0.1]}, ...
    %       'real_doppler_profile', [0,1000,0.94], ...
    %       'is_refractive_effects_removed', false, ...
    %       'is_generate_random_initial_estimates', true ...
    %   kalman_pll_config = get_kalman_pll_config(config);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Cache file path
    cache_file = 'kalman_pll_cache.mat';

    % Validate required fields in config
    % TODO: Work on error handling inside this validation function
    validate_config(config);

    % Validate if the sampling_interval from both
    % `discrete_wiener_model_config` and
    % `scintillation_training_data_config` cell arrays are the same. If
    % not, raises an error.
    sampling_interval = validate_sampling_interval(config.discrete_wiener_model_config{1,3}, ...
        config.scintillation_training_data_config{1,4});

    % Load or initialize kalman_pll_config
    [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file);

    kalman_pll_config = handle_kalman_pll_config(config, cache_file, sampling_interval, kalman_pll_config, is_cache_used);

    initial_estimates = handle_initial_estimates(config, kalman_pll_config);
end

function kalman_pll_config = handle_kalman_pll_config(config, cache_file, sampling_interval, kalman_pll_config, is_cache_used)
    % handle_kalman_pll_config
    % Computes or retrieves Kalman filter-based Phase-Locked Loop (PLL) settings.
    %
    % Syntax:
    %   kalman_pll_config = handle_kalman_pll_config(config, cache_file, sampling_interval, kalman_pll_config, is_cache_used)
    %
    % Description:
    %   This function either computes or retrieves from cache the Kalman filter settings required for PLL operation.
    %   It supports caching to reuse previously computed parameters for efficient performance. If caching is disabled,
    %   the function computes new settings, updates the Kalman PLL configuration structure, and saves it to the cache.
    %
    % Inputs:
    %   config - Struct containing all configuration details with the following fields:
    %       - training_scint_model: Specifies the scintillation model ('CSM', 'MFPSM', or 'none').
    %       - discrete_wiener_model_config: Cell array for LOS dynamics parameters to 
    %                                       be used by `get_discrete_wiener_model`.
    %         Example: {1, 3, 0.01, [0, 0, 1], 1}, where:
    %           * 1  - Number of frequency bands (`L`),
    %           * 3  - Third-order LOS dynamics (`M`),
    %           * 0.01 - Sampling time (`sampling_time`),
    %           * [0,0,1] - Sigma array for noise levels (`sigma_array`),
    %           * 1  - Ratio of its frequency to a reference frequency (`delta_array`).
    %       - scintillation_training_data_config: Cell array for scintillation model parameters.
    %         Example: {0.8, 0.7, 600, 0.01}, where:
    %           * 0.8  - S4 index (`S4`),
    %           * 0.7  - τ₀ (`tau0`),
    %           * 600  - Total simulation time (`simulation_time`),
    %           * 0.01 - Sampling time (`sampling_time`).
    %       - var_minimum_order: Minimum VAR (Vector Autoregressive) model order.
    %       - var_maximum_order: Maximum VAR model order.
    %       - C_over_N0_array_dBHz: Array of average C/N0 values for each frequency band (in dB-Hz).
    %       - is_refractive_effects_removed: Boolean flag indicating whether refractive effects are removed.
    %   cache_file - String specifying the file path for caching Kalman PLL settings.
    %   sampling_interval - Numeric value specifying the sampling interval (in seconds).
    %   kalman_pll_config - Struct to hold or update Kalman PLL settings.
    %   is_cache_used - Boolean flag indicating whether cached settings should be used.
    %
    % Outputs:
    %   kalman_pll_config - Struct containing the computed or retrieved Kalman filter settings, with the following fields:
    %       * F  - Full state transition matrix.
    %       * Q  - Full process noise covariance matrix.
    %       * H  - Measurement matrix.
    %       * R  - Measurement noise covariance matrix.
    %       * F_los, Q_los - LOS dynamics matrices.
    %       * F_var, Q_var - VAR model matrices (if applicable).
    %       * intercept_vector - VAR model intercept vector (if applicable).
    %       * var_states_amount - Number of VAR model states (if applicable).
    %       * var_model_order - Order of the VAR model (if applicable).
    %
    % Notes:
    %   - Caching is used to improve performance by avoiding repeated computations.
    %   - If `training_scint_model` is set to 'none', the function only computes LOS dynamics settings.
    %
    % Examples:
    %   % Example configuration for Kalman PLL:
    %   config = struct( ...
    %       'training_scint_model', 'CSM', ...
    %       'discrete_wiener_model_config', {1, 3, 0.01, [0, 0, 1], 1}, ...
    %       'scintillation_training_data_config', {0.8, 0.7, 300, 0.01}, ...
    %       'var_minimum_order', 1, ...
    %       'var_maximum_order', 6, ...
    %       'C_over_N0_array_dBHz', 35, ...
    %       'is_refractive_effects_removed', false ...
    %   );
    %   cache_file = 'kalman_pll_cache.mat';
    %   sampling_interval = 0.01;
    %   kalman_pll_config = struct();
    %   is_cache_used = false;
    %   kalman_pll_config = handle_kalman_pll_config(config, cache_file, sampling_interval, kalman_pll_config, is_cache_used);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Handle computation or retrieval of Kalman PLL configuration
    if is_cache_used
        fprintf('Using cached Kalman filter based PLL settings for %s.\n', config.training_scint_model);
    else
        fprintf('Computing Kalman filter based PLL settings for %s.\n', config.training_scint_model);
        
        % Compute the LOS dynamics model
        [F_los, Q_los] = get_discrete_wiener_model(config.discrete_wiener_model_config{:});

        if strcmp(config.training_scint_model, 'none')
            % Get the settings for the case where the autoregressive model
            % is not used. Thus, only the LOS dynamics state-space model is
            % outputed.
            kalman_pll_config.(config.training_scint_model) = struct('F_los', F_los, ...
                'Q_los', Q_los, ...
                'F_var', [], ...
                'Q_var', [], ...
                'F', F_los, ...
                'Q', Q_los, ...
                'H', [1, zeros(1, size(F_los, 1) - 1)], ...
                'R', diag(compute_phase_variances(config.C_over_N0_array_dBHz, sampling_interval)), ...
                'intercept_vector', [], ...
                'var_model_order', [], ...
                'var_states_amount', []);
        else
            % Get the settings for the case when the autoregressive model
            % is considered.
            kalman_pll_config = compute_settings( ...
                kalman_pll_config, ...
                config.training_scint_model, ...
                config.scintillation_training_data_config, ...
                config.var_minimum_order, ...
                config.var_maximum_order, ...
                config.C_over_N0_array_dBHz, ...
                sampling_interval, ...
                F_los, ...
                Q_los, ...
                config.is_refractive_effects_removed);
        end
        % Save updated kalman_pll_config to cache
        save(cache_file, 'kalman_pll_config');
    end
end

function initial_estimates = handle_initial_estimates(config, kalman_pll_config)
    % handle_initial_estimates
    % Generates initial state estimates and covariance matrices for Kalman filter-based PLL.
    %
    % Syntax:
    %   initial_estimates = handle_initial_estimates(config, kalman_pll_config)
    %
    % Description:
    %   This function generates initial estimates for the state vector and covariance matrix
    %   used by the Kalman filter. It supports two modes: 
    %   1. Generating random initial estimates based on uniform distributions.
    %   2. Using perfect estimates for the initial state.
    %   The function accounts for line-of-sight (LOS) dynamics and VAR (Vector Autoregressive) states,
    %   with appropriate variances based on the configuration.
    %
    % Inputs:
    %   config - Struct containing configuration details with the following fields:
    %       - is_generate_random_initial_estimates: Boolean flag indicating whether to generate random
    %         initial estimates based on `initial_states_distributions_boundaries`.
    %       - initial_states_distributions_boundaries: Cell array of boundaries for uniform distributions 
    %         to generate random initial estimates of the Doppler profile. Example:
    %           {{[-pi, pi]}, {[-5, 5]}, {[-0.1, 0.1]}, {[-0.001, 0.001]}}
    %       - real_doppler_profile: Vector containing the real Doppler profile values.
    %       - training_scint_model: String specifying the scintillation model ('CSM', 'MFPSM', or 'none').
    %
    %   kalman_pll_config - Struct containing Kalman filter settings, which must include:
    %       - var_states_amount: Number of VAR model states.
    %       - var_model_order: Order of the VAR model.
    %
    % Outputs:
    %   initial_estimates - Struct containing:
    %       * x_hat_init: Initial state vector estimate, with LOS dynamics states followed by VAR states.
    %       * P_hat_init: Initial covariance matrix, with LOS variances and uniform variance for VAR states.
    %
    % Notes:
    %   - LOS dynamics variances are computed using the uniform distributions specified in
    %     `initial_states_distributions_boundaries`.
    %   - The variance for VAR phase states is assumed to be \((\pi^2 / 3)\), corresponding to a uniform
    %     distribution within \([- \pi, \pi]\).
    %
    % Examples:
    %   % Example configuration for initial estimates:
    %   config = struct( ...
    %       'is_generate_random_initial_estimates', true, ...
    %       'initial_states_distributions_boundaries', {{[-pi, pi], [-5, 5], [-0.1, 0.1]}}, ...
    %       'real_doppler_profile', [0.01, -0.02, 0.03], ...
    %       'training_scint_model', 'CSM' ...
    %   );
    %   kalman_pll_config = struct( ...
    %       'CSM', struct('var_states_amount', 2, 'var_model_order', 1) ...
    %   );
    %   initial_estimates = handle_initial_estimates(config, kalman_pll_config);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com
    initial_estimates = struct('x_hat_init', [], 'P_hat_init', []);
    if config.is_generate_random_initial_estimates
        %%%% Generate an error vector for the initial state estimates.
        random_doppler_profile_error = arrayfun(@(col) ...
            unifrnd(config.initial_states_distributions_boundaries{1, col}(1), ...
                    config.initial_states_distributions_boundaries{1, col}(2)), ...
            1:length(config.real_doppler_profile));

        %%% Create an initial estimates vector where the first elements
        % corresponds to the line-of-sight dynamics states and the other
        % ones to the VAR states. The VAR initial states are assumed to be
        % zero arbitrarily. Nevertheless, the initial values for the VAR
        % phase estimates are not very much important, since they can be
        % any number within the range of [-pi,pi]. If we want to refactor
        % this code to be used by an extended Kalman filter instead, we
        % need to refactor this part to comprehend the amplitude estimates.
        initial_estimates.x_hat_init = [ ...
            config.real_doppler_profile.' - random_doppler_profile_error.'; ...
            zeros(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order, 1) ...
        ];

        %%% Compute line-of-sight dynamics variances according to the
        % distributions from
        % `config.initial_states_distributions_boundaries`.
        los_variances = ...
        cellfun( ...
            @(nested) ...
            (nested(2) - nested(1))^2 / 12, ...
            config.initial_states_distributions_boundaries ...
        );
        
        %%% Construct diagonal covariance matrix
        % All VAR phase states are assumed to present a variance of (pi^2/3),
        % which refers to the variance of a uniform distribution that
        % lies within the range [-pi,pi]. 
        initial_estimates.P_hat_init = ...
        blkdiag( ...
            diag(los_variances), ...
            (pi^2/3) * eye(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order) ...
        );
    else
        % Build x_hat_init with the perfect estimates for the
        % line-of-sight dynamics
        initial_estimates.x_hat_init = [ ...
            config.real_doppler_profile.'; ...
            zeros(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order, 1)
        ];
        % Build P_hat_init considering that the doppler profile is
        % perfectly known, i.e., they have covariances elements equal to 
        % zero, with the exception of the actual line-of-sight and VAR phases.
        initial_estimates.P_hat_init = blkdiag( ...
            diag([(pi^2/3),zeros(1,length(config.real_doppler_profile) - 1)]), ...
            (pi^2/3) * eye(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order) ...
        );
    end
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
    %
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
        if isfield(kalman_pll_config, config.training_scint_model) && ...
           ~isempty(fieldnames(kalman_pll_config.(config.training_scint_model)))
            if config.is_use_cached_settings
                fprintf('Using cached values for %s.\n', config.training_scint_model);
                is_cache_used = true; % Indicate cache was used
                return; % Use cached values
            else
                fprintf('Recomputing values for %s and updating cache.\n', config.training_scint_model);
            end
        else
            fprintf(['Cache found but missing settings inside %s. ' ...
                'Computing and caching new values.\n'], config.training_scint_model);
        end
    else
        fprintf('No cache file found. Initializing kalman_pll_config.\n');
        % Initialize the kalman_pll_config struct
        kalman_pll_config = struct('CSM', struct(), 'MFPSM', struct(), 'none', struct());
        bp = 1;
    end
end

function kalman_pll_config = compute_settings(kalman_pll_config, ...
    training_scint_model, scintillation_training_data_config, var_minimum_order, ...
    var_maximum_order, C_over_N0_array_dBHz, sampling_interval, ...
    F_los, Q_los, is_refractive_effects_removed)
    % compute_settings
    % Computes the Kalman filter settings and VAR model matrices based on
    % the selected scintillation model and input configuration.
    %
    % Syntax:
    %   [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = ...
    %       compute_settings(training_scint_model, scintillation_training_data_config, var_minimum_order, ...
    %       var_maximum_order, C_over_N0_array_dBHz, sampling_interval, F_los, Q_los, is_refractive_effects_removed)
    %
    % Description:
    %   This function calculates the necessary state-space matrices for the Kalman filter and
    %   VAR model parameters, including state transition matrices (F_var, F), process noise
    %   covariance matrices (Q_var, Q), measurement matrices (H), and measurement noise 
    %   covariance matrices (R). The intercept vector and the structure of the VAR model 
    %   (number of states and order) are also computed.
    %
    % Inputs:
    %   training_scint_model                      - Selected scintillation model ('CSM', 'MFPSM').
    %   scintillation_training_data_config - Cell array of scintillation model settings, such as:
    %                                         {S4, tau0, simulation_time, sampling_interval}.
    %   var_minimum_order                - Minimum VAR model order.
    %   var_maximum_order                - Maximum VAR model order.
    %   C_over_N0_array_dBHz             - Array of carrier-to-noise density ratios in dB-Hz.
    %   sampling_interval                - Sampling interval for the Kalman filter in seconds.
    %   F_los, Q_los                     - LOS dynamics state transition and covariance matrices.
    %   is_refractive_effects_removed    - Boolean flag indicating whether to remove refractive effects for MFPSM.
    %
    % Outputs:
    %   kalman_pll_config - Struct containing the computed Kalman filter 
    %                       settings, with the following fields:
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
    %   - The function preprocesses training data based on the selected scintillation model
    %     and removes refractive effects if specified (for MFPSM only).
    %   - VAR model parameters are computed using the `arfit` function.
    %   - Full Kalman matrices are constructed by combining LOS and VAR model dynamics.
    %
    % Examples:
    %   % Compute settings for the 'CSM' model:
    %   [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = ...
    %       compute_settings('CSM', {0.8, 0.7, 300, 0.01}, 1, 6, [35], 0.01, F_los, Q_los, false);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Preprocess the training data based on the scintillation model
    training_data = preprocess_training_data(training_scint_model, ...
        scintillation_training_data_config, is_refractive_effects_removed);

    % Fit VAR model
    [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = ...
        arfit(training_data, var_minimum_order, var_maximum_order);

    % Construct state transition and process noise covariance matrices for
    % the fitted VAR model
    [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices( ...
        var_coefficient_matrices, var_covariance_matrices);

    % Construct full Kalman filter matrices
    [F, Q, H, R] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, ...
        var_states_amount, var_model_order, C_over_N0_array_dBHz, sampling_interval);
    % Store results in `kalman_pll_config` struct
    kalman_pll_config.(training_scint_model) = struct('F_los', F_los, ...
        'Q_los', Q_los, ...
        'F_var', F_var, ...
        'Q_var', Q_var, ...
        'F', F, ...
        'Q', Q, ...
        'H', H, ...
        'R', R, ...
        'intercept_vector', intercept_vector, ...
        'var_model_order', var_model_order, ...
        'var_states_amount', var_states_amount);
end

function training_data = preprocess_training_data(training_scint_model, scintillation_training_data_config, is_refractive_effects_removed)
    % preprocess_training_data
    % Prepares the training data based on the selected scintillation model.
    %
    % Syntax:
    %   training_data = preprocess_training_data(training_scint_model, scintillation_training_data_config, is_refractive_effects_removed)
    %
    % Description:
    %   This function generates training data for a given scintillation model. It extracts the 
    %   phase information from the scintillation field, optionally removing refractive effects 
    %   for the MFPSM model.
    %
    % Inputs:
    %   training_scint_model                      - Scintillation model ('CSM', 'MFPSM').
    %   scintillation_training_data_config - Cell array containing scintillation model settings.
    %   is_refractive_effects_removed    - Boolean indicating whether to remove refractive effects (applies to MFPSM only).
    %
    % Outputs:
    %   training_data - Phase data extracted from the scintillation field.
    %
    % Notes:
    %   - The function calls external functions (`get_csm_data` or `get_mfpsm_data`) to generate 
    %     scintillation fields.
    %   - For MFPSM, refractive effects are removed by multiplying with the complex conjugate 
    %     of the phase screen realization.
    %
    % Examples:
    %   % Preprocess training data for the CSM model:
    %   training_data = preprocess_training_data('CSM', {0.8, 0.7, 600, 0.01}, false);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    switch training_scint_model
        case 'CSM'
            scint_complex_field = get_csm_data(scintillation_training_data_config{:});
            training_data = angle(scint_complex_field); % Extract phase data

        case 'MFPSM'
            [scint_complex_field, ps_realization] = ...
                get_mfpsm_data(scintillation_training_data_config{:});
            if is_refractive_effects_removed
                % Remove refractive effects
                scint_complex_field = scint_complex_field .* exp(-1j * ps_realization);
            end
            training_data = angle(scint_complex_field); % Extract phase data

        otherwise
            error('Unsupported scintillation model: %s', training_scint_model);
    end
end

function [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(var_coefficient_matrices, var_covariance_matrices)
    % construct_var_matrices
    % Constructs VAR model matrices from coefficient and covariance matrices.
    %
    % Syntax:
    %   [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(var_coefficient_matrices, var_covariance_matrices)
    %
    % Description:
    %   This function generates the state transition matrix (F_var) and the process 
    %   noise covariance matrix (Q_var) for a VAR (Vector Autoregressive) model.
    %
    % Inputs:
    %   var_coefficient_matrices - Coefficient matrices for the VAR model.
    %   var_covariance_matrices  - Covariance matrices for the VAR model.
    %
    % Outputs:
    %   F_var          - State transition matrix for the VAR model.
    %   Q_var          - Process noise covariance matrix for the VAR model.
    %   var_states_amount - Number of states in the VAR model.
    %   var_model_order  - Order of the VAR model.
    %
    % Notes:
    %   - The function constructs augmented matrices to handle higher-order VAR models.
    %   - Q_var includes zero padding for augmented states.
    %
    % Examples:
    %   [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(coeff_matrix, cov_matrix);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    var_states_amount = size(var_coefficient_matrices, 1);
    var_model_order = size(var_coefficient_matrices, 2) / var_states_amount;

    % Construct augmented state transition matrix
    F_var = [var_coefficient_matrices; [eye(var_states_amount * ...
        (var_model_order - 1)), zeros(var_states_amount * ...
        (var_model_order - 1), var_states_amount)]];

    % Construct augmented process noise covariance matrix
    Q_var = blkdiag(var_covariance_matrices, zeros(var_states_amount * ...
        (var_model_order - 1)));
end

function [F, Q, H, R] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, var_states_amount, var_model_order, C_over_N0_array_dBHz, sampling_interval)
    % construct_kalman_matrices
    % Constructs the full Kalman filter matrices.
    %
    % Syntax:
    %   [F, Q, H, R] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, var_states_amount, var_model_order, C_over_N0_array_dBHz, sampling_interval)
    %
    % Description:
    %   This function combines LOS dynamics and VAR model matrices to construct 
    %   the state transition matrix (F), process noise covariance matrix (Q), 
    %   measurement matrix (H), and measurement noise covariance matrix (R) 
    %   for the Kalman filter.
    %
    % Inputs:
    %   F_los, Q_los     - LOS dynamics state transition and covariance matrices.
    %   F_var, Q_var     - VAR model state transition and covariance matrices.
    %   var_states_amount - Number of states in the VAR model.
    %   var_model_order  - Order of the VAR model.
    %   C_over_N0_array_dBHz - Array of C/N0 values in dB-Hz.
    %   sampling_interval - Sampling interval for computations.
    %
    % Outputs:
    %   F - Full state transition matrix.
    %   Q - Full process noise covariance matrix.
    %   H - Measurement matrix.
    %   R - Measurement noise covariance matrix.
    %
    % Notes:
    %   - The measurement matrix (H) maps the LOS and VAR states to measurements.
    %   - Measurement noise covariance (R) is computed from the provided C/N0 values.
    %
    % Examples:
    %   [F, Q, H, R] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, 3, 2, [35], 0.01);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    F = blkdiag(F_los, F_var); % Combine LOS and VAR state transitions
    Q = blkdiag(Q_los, Q_var); % Combine LOS and VAR covariance matrices
    H = [1, zeros(1, size(F_los, 1) - 1), 1, ... 
         zeros(1, var_states_amount * var_model_order - 1)]; % Measurement matrix
    R = diag(compute_phase_variances(C_over_N0_array_dBHz, sampling_interval)); % Measurement noise covariance
end

function sigma2_array = compute_phase_variances(C_over_N0_array_dBHz, sampling_interval)
    % compute_phase_variances
    % Computes the measurement noise variances based on C/N0 values.
    %
    % Syntax:
    %   sigma2_array = compute_phase_variances(C_over_N0_array_dBHz, sampling_interval)
    %
    % Description:
    %   This function calculates the phase noise variances (sigma^2) in radians squared 
    %   based on the carrier-to-noise ratio (C/N0) in dB-Hz and the sampling interval in seconds.
    %
    % Inputs:
    %   C_over_N0_array_dBHz - Array of carrier-to-noise density ratios in dB-Hz.
    %   sampling_interval    - Sampling interval in seconds.
    %
    % Outputs:
    %   sigma2_array - Array of computed measurement noise variances.
    %
    % Notes:
    %   - The calculation uses standard relationships between C/N0 and phase noise.
    %
    % Examples:
    %   sigma2_array = compute_phase_variances([35, 40], 0.01);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com
    
    C_over_N0_array_linear = 10.^(C_over_N0_array_dBHz ./ 10); % Convert dB-Hz to linear scale
    sigma2_array = (1 ./ (2 * C_over_N0_array_linear * sampling_interval)) ...
        .* (1 + 1 ./ (2 * C_over_N0_array_linear * sampling_interval)); % Compute variances
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
    %
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

function validate_config(config)
    % TODO: Work on error handling inside this validation function
    % Ensure required fields are present in the configuration struct
    required_fields = {'discrete_wiener_model_config', ...
                       'scintillation_training_data_config', ...
                       'var_minimum_order', ...
                       'var_maximum_order', ...
                       'C_over_N0_array_dBHz', ...
                       'training_scint_model', ...
                       'initial_states_distributions_boundaries', ...
                       'real_doppler_profile', ...
                       'is_refractive_effects_removed', ...
                       'is_use_cached_settings', ...
                       'is_generate_random_initial_estimates' ...
                       };
    missing_fields = setdiff(required_fields, fieldnames(config));
    if ~isempty(missing_fields)
        error('get_kalman_pll_config:MissingFields', ...
              'The following required fields are missing: %s', strjoin(missing_fields, ', '));
    end
end

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
    %       - training_scint_model: Specifies the scintillation model ('CSM', 'TPPSM', or 'none').
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
                'W', zeros(size(F_los,1),1), ...
                'var_model_order', [], ...
                'var_states_amount', []);
        else
            % Get the settings for the case when the autoregressive model
            % is considered.
            kalman_pll_config = compute_settings( ...
                kalman_pll_config, ...
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
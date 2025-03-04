function kalman_pll_config = update_cache(general_config, cache_file, kalman_pll_config, is_cache_used, is_enable_cmd_print)
% update_cache
%
% Syntax:
%   kalman_pll_config = update_cache(general_config, cache_file, kalman_pll_config, is_cache_used)
%
% Description:
%   Either computes new Kalman PLL settings or retrieves them from cache for efficiency.
%   If caching is disabled, the function computes new settings, updates the Kalman PLL configuration
%   struct, and saves it to the specified cache file.
%
% Inputs:
%   general_config - Struct containing all configuration details with required fields:
%       - discrete_wiener_model_config: Cell array for LOS dynamics parameters to be used by
%         the get_los_phase function. Example: {1, 3, 0.01, [0,0,1], 1}.
%
%       - scintillation_training_data_config: Struct with fields:
%             S4, tau0, simulation_time, and sampling_interval.
%
%       - var_minimum_order: Minimum VAR model order.
%       - var_maximum_order: Maximum VAR model order.
%       - C_over_N0_array_dBHz: Numeric vector (positive) of average C/N0 values.
%       - initial_states_distributions_boundaries: Non-empty cell array where each element is a
%         1×2 numeric vector (with the first element less than the second).
%       - real_doppler_profile: Non-empty numeric vector.
%       - augmentation_model_initializer: 
%           * 'arfit' (Initializes multivariate autoregressive model parameters. It also calculates an intercept vector.), 
%           * 'aryule' (Initializes the autoregressive model using the Yule-Walker method, OBS: Signal Processing Toolbox needed), 
%           * 'rbf', (Initializes the Radial Basis Function Network Weights); 
%       - is_use_cached_settings: Boolean flag indicating whether cached configurations should be used.
%       - is_generate_random_initial_estimates: Boolean flag indicating whether initial estimates are randomly generated.
%
%   cache_file        - String specifying the file path for caching results.
%   kalman_pll_config - Struct to hold or update the Kalman PLL settings.
%   is_cache_used     - Boolean indicating whether cached settings should be used.
%   is_enable_cmd_print - Boolean flag for enabling the command prints.
%                         It is recommended to disable this option for
%                         monte carlo runs.
% Outputs:
%   kalman_pll_config - Struct containing the computed Kalman filter settings.
%       A field is added to the struct with the name of the scintillation model (e.g., 'CSM'
%       or 'TPPSM'). This substruct has the following fields:
%           * F_los             : LOS dynamics state transition matrix.
%           * Q_los             : LOS dynamics process noise covariance matrix.
%           * F_var             : VAR model state transition matrix.
%           * Q_var             : VAR model process noise covariance matrix.
%           * F                 : Full state transition matrix.
%           * Q                 : Full process noise covariance matrix.
%           * H                 : Measurement matrix.
%           * R                 : Measurement noise covariance matrix.
%           * W                 : Additional matrix from construct_kalman_matrices.
%           * intercept_vector  : VAR model intercept vector.
%           * var_model_order   : Order of the VAR model.
%           * var_states_amount : Number of VAR model states.
%
% Example:
%   general_config = struct( ...
%       'discrete_wiener_model_config', {discrete_wiener_model_config}, ...
%       'scintillation_training_data_config', {scintillation_training_data_config}, ...
%       'arfit_model_order', var_minimum_order,
%       'C_over_N0_array_dBHz', C_over_N0_array_dBHz, ...
%       'initial_states_distributions_boundaries', {initial_states_distributions_boundaries}, ...
%       'real_doppler_profile', real_doppler_profile, ...
%       'augmentation_model_initializer', 'aryule', ...
%       'is_use_cached_settings', is_use_cached_settings, ...
%       'is_generate_random_initial_estimates', is_generate_random_initial_estimates ...
%   );
%   cache_file = 'kalman_pll_cache.mat';
%   kalman_pll_config = struct();
%   is_cache_used = false;
%   kalman_pll_config = update_cache(general_config, cache_file, kalman_pll_config, is_cache_used);
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Validate inputs
    validateattributes(kalman_pll_config, {'struct'}, {'nonempty'}, mfilename, 'kalman_pll_config');
    validateattributes(is_cache_used, {'logical'}, {'scalar'}, mfilename, 'is_cache_used');

    if is_cache_used
        if is_enable_cmd_print
            fprintf('Using cached Kalman filter-based PLL settings.\n');
        end
    else
        if is_enable_cmd_print
            fprintf('Computing Kalman filter-based PLL settings.\n');
        end

        % Compute the complete Kalman settings (VAR model, etc.)
        kalman_pll_config = build_kalman_pll_config( ...
            general_config, ...
            kalman_pll_config ...
            );

        % Save updated kalman_pll_config to the cache file
        if is_enable_cmd_print
            fprintf('Caching Kalman filter-based PLL settings.\n');
        end
        save(cache_file, 'kalman_pll_config');
    end
end

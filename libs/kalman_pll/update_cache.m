function kalman_pll_config = update_cache(general_config, cache_file, kalman_pll_config, is_cache_used)
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
%   general_config - Struct containing configuration details with required fields:
%       - discrete_wiener_model_config: Cell array for LOS dynamics parameters used by
%         get_discrete_wiener_model.
%       - scintillation_training_data_config: Struct containing scintillation model settings.
%           For CSM, expected fields:
%               scintillation_model - Must be 'CSM'
%               S4                  - Scintillation index (0 <= S4 <= 1)
%               tau0                - Signal decorrelation time (positive scalar)
%               simulation_time     - Duration of simulation (positive scalar)
%               sampling_interval   - Sampling interval (positive scalar)
%           For TPPSM, expected fields:
%               scintillation_model - Must be 'TPPSM'
%               scenario            - A string specifying the scenario ('Weak', 'Moderate', 'Severe')
%               simulation_time     - Duration of simulation (positive scalar)
%               sampling_interval   - Sampling interval (positive scalar)
%               is_refractive_effects_removed - Boolean flag (true or false)
%       - var_minimum_order: Minimum VAR model order.
%       - var_maximum_order: Maximum VAR model order.
%       - C_over_N0_array_dBHz: Array of average C/N0 values for each frequency band (in dB-Hz).
%       - initial_states_distributions_boundaries: Cell array for initial state distribution boundaries.
%       - real_doppler_profile: Array or structure with the real Doppler profile.
%       - is_use_cached_settings: Boolean flag indicating whether to use cached settings.
%       - is_generate_random_initial_estimates: Boolean flag for generating random initial estimates.
%
%   cache_file        - String specifying the file path for caching results.
%   kalman_pll_config - Struct to hold or update the Kalman PLL settings.
%   is_cache_used     - Boolean indicating whether cached settings should be used.
%
% Outputs:
%   kalman_pll_config - Struct containing the computed or retrieved Kalman filter settings
%                       (F, Q, H, R, F_los, Q_los, F_var, Q_var, etc.).
%
% Example:
%   general_config = struct( ...
%       'discrete_wiener_model_config', {discrete_wiener_model_config}, ...
%       'scintillation_training_data_config', {scintillation_training_data_config}, ...
%       'var_minimum_order', var_minimum_order, ...
%       'var_maximum_order', var_maximum_order, ...
%       'C_over_N0_array_dBHz', C_over_N0_array_dBHz, ...
%       'initial_states_distributions_boundaries', {initial_states_distributions_boundaries}, ...
%       'real_doppler_profile', real_doppler_profile, ...
%       'is_use_cached_settings', is_use_cached_settings, ...
%       'is_generate_random_initial_estimates', is_generate_random_initial_estimates ...
%   );
%   cache_file = 'kalman_pll_cache.mat';
%   kalman_pll_config = struct();
%   is_cache_used = false;
%   kalman_pll_config = update_cache(general_config, cache_file, kalman_pll_config, is_cache_used);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Validate inputs
    validateattributes(general_config, {'struct'}, {'nonempty'}, mfilename, 'general_config');
    req_config_fields = {
        'discrete_wiener_model_config', ...
        'scintillation_training_data_config', ...
        'var_minimum_order', ...
        'var_maximum_order', ...
        'C_over_N0_array_dBHz'
    };
    for fc = 1:numel(req_config_fields)
        if ~isfield(general_config, req_config_fields{fc})
            error('update_cache:MissingField', ...
                'Config struct is missing the field "%s".', req_config_fields{fc});
        end
    end
    
    validateattributes(cache_file, {'char','string'}, {'nonempty'}, mfilename, 'cache_file');
    validateattributes(kalman_pll_config, {'struct'}, {'nonempty'}, mfilename, 'kalman_pll_config');
    validateattributes(is_cache_used, {'logical'}, {'scalar'}, mfilename, 'is_cache_used');

    scint_training_data_cfg = general_config.scintillation_training_data_config;
    validateattributes(scint_training_data_cfg, {'struct'}, {'nonempty'}, mfilename, 'scintillation_training_data_config');
    if ~isfield(scint_training_data_cfg, 'sampling_interval')
        error('update_cache:MissingField', 'scintillation_training_data_config is missing the field "sampling_interval".');
    end
    % Extract and validate sampling_interval from st_config.
    sampling_interval = scint_training_data_cfg.sampling_interval;
    validateattributes(sampling_interval, {'numeric'}, {'scalar','real','positive'}, mfilename, 'sampling_interval');

    if is_cache_used
        fprintf('Using cached Kalman filter-based PLL settings.\n');
    else
        fprintf('Computing Kalman filter-based PLL settings.\n');

        % Compute the LOS dynamics model
        [F_los, Q_los] = get_discrete_wiener_model(general_config.discrete_wiener_model_config{:});

        % Compute the complete Kalman settings (VAR model, etc.)
        kalman_pll_config = build_kalman_pll_config( ...
            kalman_pll_config, ...
            scint_training_data_cfg, ...
            general_config.var_minimum_order, ...
            general_config.var_maximum_order, ...
            general_config.C_over_N0_array_dBHz, ...
            F_los, ...
            Q_los ...
            );

        % Save updated kalman_pll_config to the cache file
        fprintf('Caching Kalman filter-based PLL settings.\n');
        save(cache_file, 'kalman_pll_config');
    end
end

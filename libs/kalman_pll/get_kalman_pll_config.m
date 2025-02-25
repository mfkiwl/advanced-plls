function [kalman_pll_config, initial_estimates] = get_kalman_pll_config(general_config, cache_dir)
    % get_kalman_pll_config
    % Generates or retrieves Kalman filter settings and initial estimates,
    % based on the provided configuration parameters.
    %
    % Syntax:
    %   [kalman_pll_config, initial_estimates] = get_kalman_pll_config(general_config)
    %
    % Description:
    %   This function computes or retrieves the Kalman filter settings—including state‐space
    %   matrices, noise covariance matrices, and initial estimates—based on the provided configuration.
    %   It uses caching to avoid repeated computations. If caching is disabled or the cache is invalid,
    %   new parameters are computed and the cache is updated.
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
    %       - is_use_cached_settings: Boolean flag indicating whether cached configurations should be used.
    %       - is_generate_random_initial_estimates: Boolean flag indicating whether initial estimates are randomly generated.
    %   cache_dir - directory of the cache folder
    %
    % Outputs:
    %   kalman_pll_config - Struct containing the computed Kalman filter settings, with fields:
    %           F, Q, H, R, F_los, Q_los, F_var, Q_var, intercept_vector,
    %           var_states_amount, and var_model_order.
    %   initial_estimates - Column vector of initial estimates. Its number of rows equals the number of states in the state transition matrix.
    %
    % Notes:
    %   - The sampling interval is taken from scintillation_training_data_config.sampling_interval,
    %     and it is compared to the sampling interval in discrete_wiener_model_config (element 3).
    %   - Lower-level functions (get_cached_kalman_pll_config, update_cache, and get_initial_estimates)
    %     handle further validations.
    %
    % Examples:
    %   general_config = struct( ...
    %       'discrete_wiener_model_config', { {1,3,0.01,[0,0,1e-2],1} }, ...
    %       'scintillation_training_data_config', struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0',0.7, 'simulation_time',300, 'sampling_interval',0.01), ...
    %       'var_minimum_order', 1, ...
    %       'var_maximum_order', 6, ...
    %       'C_over_N0_array_dBHz', [35], ...
    %       'initial_states_distributions_boundaries',{ {[-pi,pi], [-5,5], [-0.1,0.1]} }, ...
    %       'real_doppler_profile', [0,1000,0.94], ...
    %       'is_use_cached_settings', false, ...
    %       'is_generate_random_initial_estimates', true ...
    %   );
    %   [kalman_pll_config, initial_estimates] = get_kalman_pll_config(general_config);
    %
    % Author: Rodrigo de Lima Florindo
    % ORCID: https://orcid.org/0000-0003-0412-5583
    % Email: rdlfresearch@gmail.com

    % Validate required fields in general_config
    required_fields = {'discrete_wiener_model_config', 'scintillation_training_data_config', ...
                       'var_minimum_order', 'var_maximum_order', 'C_over_N0_array_dBHz', ...
                       'initial_states_distributions_boundaries', 'real_doppler_profile', ...
                       'is_use_cached_settings', 'is_generate_random_initial_estimates'};
    for i = 1:numel(required_fields)
        if ~isfield(general_config, required_fields{i})
            error('get_kalman_pll_config:MissingField', ...
                'general_config is missing the required field: %s', required_fields{i});
        end
    end

    % Validate sampling_interval in scintillation_training_data_config
    st_config = general_config.scintillation_training_data_config;
    if ~isfield(st_config, 'sampling_interval') || isempty(st_config.sampling_interval)
        error('get_kalman_pll_config:MissingField', ...
            'scintillation_training_data_config must have a non-empty sampling_interval field.');
    end
    sampling_interval = st_config.sampling_interval;
    
    % Validate initial_states_distributions_boundaries: must be a non-empty cell array of 1x2 numeric vectors with min < max.
    if ~iscell(general_config.initial_states_distributions_boundaries) || isempty(general_config.initial_states_distributions_boundaries)
        error('get_kalman_pll_config:InvalidBoundaries', 'initial_states_distributions_boundaries must be a non-empty cell array.');
    end
    for i = 1:length(general_config.initial_states_distributions_boundaries)
        b = general_config.initial_states_distributions_boundaries{i};
        validateattributes(b, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'initial_states_distributions_boundaries');
        if b(1) >= b(2)
            error('get_kalman_pll_config:InvalidBoundaries', 'Each boundary must have its first element less than its second.');
        end
    end

    % Validate real_doppler_profile: non-empty numeric vector.
    validateattributes(general_config.real_doppler_profile, {'numeric'}, {'nonempty','vector'}, mfilename, 'real_doppler_profile');

    % Validate C_over_N0_array_dBHz: non-empty, positive numeric vector.
    validateattributes(general_config.C_over_N0_array_dBHz, {'numeric'}, {'nonempty','vector','positive'}, mfilename, 'C_over_N0_array_dBHz');
    
    % Validate `cache_file` path
    validateattributes(cache_dir, {'char', 'string'}, {'nonempty'}, mfilename, 'cache_file');

    % Compare sampling_interval in discrete_wiener_model_config and scintillation_training_data_config
    dw_config = general_config.discrete_wiener_model_config;
    if numel(dw_config) < 3
        error('get_kalman_pll_config:MissingField', 'discrete_wiener_model_config must have at least 3 elements.');
    end
    sampling_interval_dw = dw_config{3};
    validateattributes(sampling_interval_dw, {'numeric'}, {'scalar','real','positive'}, mfilename, 'discrete_wiener_model_config{3}');
    if abs(sampling_interval - sampling_interval_dw) > 1e-10
        error('get_kalman_pll_config:SamplingIntervalMismatch', ...
            'Sampling intervals in scintillation_training_data_config and discrete_wiener_model_config must match.');
    end

    % Ensure the cache directory exists
    if ~exist(cache_dir, 'dir')
        mkdir(cache_dir);
    end
    cache_file = fullfile(cache_dir, 'kalman_pll_cache.mat');

    [kalman_pll_config, is_cache_used] = get_cached_kalman_pll_config(general_config, cache_file);
    kalman_pll_config = update_cache(general_config, cache_file, kalman_pll_config, is_cache_used);
    initial_estimates = get_initial_estimates(general_config, kalman_pll_config);
end

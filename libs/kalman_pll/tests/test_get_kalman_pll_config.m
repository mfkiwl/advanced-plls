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
    %       - training_scint_model: Selected scintillation model ('CSM', 'TPPSM', or 'none').
    %       - initial_states_distributions_boundaries: Uniform
    %          distributions boundaries for generating the initial state
    %          estimates for the Doppler profile used by the Kalman Filter.
    %       - real_doppler_profile: Real Doppler profile used to simulate the
    %          synthetic line-of-sight dynamics.
    %       - is_refractive_effects_removed: Boolean flag indicating whether refractive 
    %          effects are removed for TPPSM.
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

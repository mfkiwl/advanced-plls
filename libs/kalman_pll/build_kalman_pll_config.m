function kalman_pll_config = build_kalman_pll_config(general_config, ...
    kalman_pll_config)
% build_kalman_pll_config
%
% Computes the Kalman filter settings and VAR model matrices based on the
% selected scintillation model and input configuration.
%
% Syntax:
%   kalman_pll_config = build_kalman_pll_config(kalman_pll_config, ...
%       scintillation_training_data_config, var_minimum_order, ...
%       var_maximum_order, C_over_N0_array_dBHz, F_los, Q_los)
%
% Description:
%   This function calculates the necessary state-space matrices for the Kalman
%   filter and VAR model parameters, including:
%       - F_var, F: State transition matrices for the VAR model and full system.
%       - Q_var, Q: Process noise covariance matrices for the VAR model and full system.
%       - H: Measurement matrix.
%       - R: Measurement noise covariance matrix.
%       - W: An additional matrix from construct_kalman_matrices.
%       - intercept_vector: VAR model intercept vector.
%       - var_states_amount: Number of VAR model states.
%       - var_model_order: Order of the VAR model.
%
%   The sampling_interval is extracted from scintillation_training_data_config.
%
% Inputs:
%   kalman_pll_config - Struct to hold or update the Kalman PLL settings.
%       It is expected to be a struct that either is empty (if no settings
%       have been computed yet) or contains previously computed settings.
%       For example, if the "CSM" scintillation model has been computed before,
%       kalman_pll_config might have a field "CSM" with subfields:
%           F_los, Q_los, F_var, Q_var, F, Q, H, R, W, intercept_vector,
%           var_model_order, and var_states_amount.
%
%   general_config - Struct containing configuration settings:
%       discrete_wiener_model_config: Cell array {L, M, sampling_interval, sigma, delta}
%           L        - Number of carriers (positive integer scalar)
%           M        - Order of the Wiener process (positive integer scalar)
%           sampling_interval - Sampling interval for LOS dynamics (positive scalar)
%           sigma    - Numeric vector of variances; length must equal (L + M - 1)
%           delta    - Numeric vector of normalized frequencies; length must equal L
%
%       scintillation_training_data_config: Struct with scintillation model settings.
%           scintillation_model - String indicating the model ('CSM', 'TPPSM', or 'NONE')
%           For CSM:
%               S4              - Scintillation index (numeric scalar in [0,1])
%               tau0            - Signal decorrelation time (positive scalar)
%               simulation_time - Duration of simulation (positive scalar)
%               sampling_interval - Sampling interval (positive scalar)
%           For TPPSM:
%               scenario        - String: 'Weak', 'Moderate', or 'Severe'
%               simulation_time - Duration of simulation (positive scalar)
%               sampling_interval - Sampling interval (positive scalar)
%               is_refractive_effects_removed - Logical flag (true/false; defaults to false)
%
%       C_over_N0_array_dBHz - Numeric vector (positive) of average C/N0 values (in dBHz)
%       initial_states_distributions_boundaries - Non-empty cell array; each cell contains a 1x2 numeric vector
%           specifying lower and upper bounds (first element < second element).
%       real_doppler_profile - Non-empty numeric vector of Doppler profile values.
%       augmentation_model_initializer - Struct specifying the augmentation model initialization method and its parameters:
%           Fields:
%               id - A string indicating the initialization method. Allowed values are: 'arfit', 'aryule', 'rbf', or 'none'.
%               model_params - A struct containing method-specific parameters:
%                   * For 'arfit' and 'aryule': must include the field 'model_order' (a numeric value).
%                   * For 'rbf': must include the field 'neurons_amount' (a numeric value).
%
%       is_use_cached_settings - Boolean flag to use cached configurations if available.
%       is_generate_random_initial_estimates - Boolean flag to generate initial estimates randomly.
%
%   F_los, Q_los          - LOS dynamics state transition and covariance matrices.
%
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
%   kalman_pll_config = build_kalman_pll_config(kalman_pll_config, ...
%       scintillation_training_data_config, 1, 6, [35], F_los, Q_los);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Compute the LOS dynamics model
    [F_los, Q_los] = get_discrete_wiener_model(general_config.discrete_wiener_model_config{:});

    % Preprocess Training Data
    training_data = preprocess_training_data(general_config.scintillation_training_data_config);

    switch general_config.augmentation_model_initializer.id
        case 'arfit'
            [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = ...
                arfit(training_data, general_config.augmentation_model_initializer.model_params.model_order, general_config.augmentation_model_initializer.model_params.model_order);
            % Construct State Transition and Process Noise Matrices
            [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices( ...
                var_coefficient_matrices, var_covariance_matrices);
        
            % Construct Full Kalman Filter Matrices
            [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, ...
                intercept_vector, var_states_amount, var_model_order, ...
                general_config.C_over_N0_array_dBHz, general_config.scintillation_training_data_config.sampling_interval);
        case 'aryule'
            [ar_coefficients, ar_variance] = aryule(training_data, general_config.augmentation_model_initializer.model_params.model_order);
            % Construct State Transition and Process Noise Matrices
            [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices( ...
                ar_coefficients, ar_variance);
        
            % Construct Full Kalman Filter Matrices
            [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, ...
                [], var_states_amount, var_model_order, ...
                general_config.C_over_N0_array_dBHz, general_config.scintillation_training_data_config.sampling_interval);
        case 'rbf'
            % NOTE: For now, does nothing. This part is under development.
            % This case well never happen, considering that there is an
            % error in get_kalman_pll_config.m that is rasing to prevent it.
        case 'none'
            % For 'none', no augmentation is applied.
            % In this simple implementation, we leave the related matrices undefined.
            var_model_order = [];
            var_states_amount = [];
            F_var = [];
            Q_var = [];
            F = F_los;
            Q = Q_los;
            H = [1, zeros(1, size(F_los,1)-1)];
            R = diag(compute_phase_variances(general_config.C_over_N0_array_dBHz, general_config.discrete_wiener_model_config{3}));
            W = zeros(size(F_los,1), 1);
    end

    % Store Results in Output Struct
    kalman_pll_config.(general_config.scintillation_training_data_config.scintillation_model) = struct( ...
        'F_los', F_los, ...
        'Q_los', Q_los, ...
        'F_var', F_var, ...
        'Q_var', Q_var, ...
        'F', F, ...
        'Q', Q, ...
        'H', H, ...
        'R', R, ...
        'W', W, ...
        'var_model_order', var_model_order, ...
        'var_states_amount', var_states_amount ...
    );
end

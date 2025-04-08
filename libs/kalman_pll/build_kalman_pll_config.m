function kalman_pll_config = build_kalman_pll_config(general_config, ...
    kalman_pll_config)
% build_kalman_pll_config
%
% Summary:
%   Configures and constructs the Kalman filter and VAR (or alternative
%   augmentation model) state-space matrices based on the provided configuration.
%   The function computes line-of-sight (LOS) dynamics using a discrete Wiener model,
%   preprocesses training data, and depending on the selected augmentation model,
%   computes the necessary matrices for the state-space filter.
%
% Syntax:
%   kalman_pll_config = build_kalman_pll_config(general_config, kalman_pll_config)
%
% Description:
%   This function derives the Kalman filter and VAR/augmentation model settings 
%   from user-specified configuration parameters. It performs the following steps:
%
%     1. Computes the LOS dynamics matrices (F_los and Q_los) by calling the discrete
%        Wiener model with parameters specified in general_config.discrete_wiener_model_config.
%     2. Preprocesses the training data extracted from general_config.scintillation_training_data_config.
%     3. Depending on the augmentation model identifier provided in
%        general_config.augmentation_model_initializer.id, it computes additional system
%        matrices:
%           - 'arfit': Estimates a VAR model using ARFIT to obtain the state transition (F_var)
%                      and process noise (Q_var) matrices, and then constructs the full Kalman
%                      filter matrices.
%           - 'aryule': Similar to 'arfit', but uses the ARYULE approach for parameter estimation.
%           - 'kinematic': Constructs an augmented model for single-frequency tracking using an
%                          alternative LOS augmentation based on a discrete Wiener model.
%           - 'arima': Constructs an ARIMA state-space model by estimating ARIMA parameters and
%                      embedding the model into the Kalman filter framework.
%           - 'rbf': Reserved for a radial basis function approach (currently under development).
%           - 'none': Uses the LOS dynamics matrices alone, without any state augmentation.
%
%   Finally, the computed matrices and additional settings are stored in the output
%   structure under a field named after the selected scintillation model (e.g., 'CSM', 'TPPSM').
%
% Inputs:
%   general_config - Struct containing configuration settings including:
%
%       discrete_wiener_model_config: Cell array {L, M, sampling_interval, sigma, delta}
%           L                 - Number of carriers (positive integer scalar).
%           M                 - Order of the Wiener process (positive integer scalar).
%           sampling_interval - Sampling interval for LOS dynamics (positive scalar).
%           sigma             - Numeric vector of variances (length must equal L + M - 1).
%           delta             - Numeric vector of normalized frequencies (length must equal L).
%
%       scintillation_training_data_config: Struct with scintillation model settings:
%           scintillation_model - String indicating the model ('CSM', 'TPPSM', or 'NONE').
%           For 'CSM':
%               S4              - Scintillation index (numeric scalar in [0,1]).
%               tau0            - Signal decorrelation time (positive scalar).
%               simulation_time - Total simulation duration (positive scalar).
%               sampling_interval - Sampling interval (positive scalar).
%           For 'TPPSM':
%               scenario        - String: 'Weak', 'Moderate', or 'Severe'.
%               simulation_time - Total simulation duration (positive scalar).
%               sampling_interval - Sampling interval (positive scalar).
%               is_refractive_effects_removed - Logical flag (true/false; defaults to false).
%
%       C_over_N0_array_dBHz - Numeric vector (positive) of average C/N0 values in dBHz.
%
%       initial_states_distributions_boundaries - Cell array where each cell contains a 1x2
%           numeric vector specifying lower and upper bounds for initial state estimates.
%
%       real_doppler_profile - Numeric vector containing Doppler profile values.
%
%       augmentation_model_initializer - Struct with fields for model initialization:
%           id         - String identifier of the augmentation model. Allowed values:
%                        'arfit', 'aryule', 'kinematic', 'arima', 'rbf', or 'none'.
%           model_params - Struct with method-specific parameters:
%                           * For 'arfit' and 'aryule': expects 'model_order' (numeric).
%                           * For 'kinematic': expects 'wiener_mdl_order' and
%                                            'process_noise_variance'.
%                           * For 'arima': expects 'p', 'D', and 'q'.
%                           * For 'rbf': expects 'neurons_amount'.
%
%       is_use_cached_settings - Boolean flag indicating if cached settings should be reused.
%       is_generate_random_initial_estimates - Boolean flag for random generation of initial estimates.
%
%   kalman_pll_config - Struct that holds or updates the Kalman PLL settings. If empty,
%                       new settings are computed; otherwise, the existing configuration
%                       is updated.
%
% Note:
%   The LOS dynamics matrices (F_los and Q_los) are internally computed from the
%   discrete_wiener_model_config using get_discrete_wiener_model, so these do not have
%   to be provided as input.
%
% Outputs:
%   kalman_pll_config - Struct that includes a substructure named after the selected
%                       scintillation model. This substructure comprises:
%           F   : Full state transition matrix.
%           Q   : Full process noise covariance matrix.
%           H   : Measurement matrix.
%           R   : Measurement noise covariance matrix.
%           W   : Additional matrix derived from construct_kalman_matrices.
%           augmentation_model_initializer:
%                The original augmentation model initializer structure.
%
% Example:
%   kalman_pll_config = build_kalman_pll_config(general_config, kalman_pll_config);
%
% References:
%   [1] Durbin, J., & Koopman, S. J. (2001). Time Series Analysis by State Space 
%       Methods.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Compute the LOS dynamics model
    [F_los, Q_los] = get_discrete_wiener_model(general_config.discrete_wiener_model_config{:});

    % Preprocess Training Data
    training_data = preprocess_training_data(general_config.scintillation_training_data_config);
    sampling_interval = general_config.discrete_wiener_model_config{3};

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
                general_config.C_over_N0_array_dBHz, sampling_interval);
        case 'aryule'
            [ar_coefficients, ar_variance] = aryule(training_data, general_config.augmentation_model_initializer.model_params.model_order);
            % Construct State Transition and Process Noise Matrices
            [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices( ...
                -ar_coefficients(2:end), ar_variance);
        
            % Construct Full Kalman Filter Matrices
            [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, ...
                0, var_states_amount, var_model_order, ...
                general_config.C_over_N0_array_dBHz, sampling_interval);
        case 'kinematic'
            L_aug = 1; % Single-frequency tracking
            M_aug = general_config.augmentation_model_initializer.model_params.wiener_mdl_order;
            sampling_interval = general_config.discrete_wiener_model_config{3};
            sigma_aug = [0,0,general_config.augmentation_model_initializer.model_params.process_noise_variance];
            delta_array_aug = 1;
            [F_wiener_aug, Q_wiener_aug] = get_discrete_wiener_model(L_aug, M_aug, sampling_interval, sigma_aug, delta_array_aug);
            F = blkdiag(F_los,F_wiener_aug);
            Q = blkdiag(Q_los,Q_wiener_aug);
            H = [1, zeros(1, size(F_los,1) - 1), 1, zeros(1, size(F_wiener_aug,1)-1)];
            R = diag(compute_phase_variances(general_config.C_over_N0_array_dBHz, general_config.discrete_wiener_model_config{3}));
            W = zeros(size(F_los,1) + size(F_wiener_aug,1), 1);
        case 'arima'
            % The development of the matrices was inspired in what was
            % proposed in [Section 3.4, 1]

            p = general_config.augmentation_model_initializer.model_params.p;
            D = general_config.augmentation_model_initializer.model_params.D;
            q = general_config.augmentation_model_initializer.model_params.q;
            % Model parameter estimation
            mdl = arima(p,D,q);
            mdl.Constant = 0;
            est_mdl = estimate(mdl,training_data);
            r = max(p,q + 1);
            
            % ARMA State Space Model (SSM):
            ar_params = cell2mat(est_mdl.AR).';
            ma_params = cell2mat(est_mdl.MA).';
            variance = est_mdl.Variance;

            phi_vec = zeros(r,1);
            phi_vec(1:length(ar_params)) = ar_params;

            theta_vec = zeros(1,r-1);
            theta_vec(1:length(ma_params)) = ma_params;
            F_arma = [phi_vec,[eye(r-1);zeros(1,r-1)]];
            G_arma = [1,theta_vec].';

            % ARIMA SSM:
            % Note: The ARIMA SSM can be seen as an extension of the ARMA
            % SSM, including also the previous values of the states differences
            upper_right_matrix = zeros(D,r);
            upper_right_matrix(:,1) = 1;
            F_arima = [triu(ones(D)), upper_right_matrix; [zeros(r,D), F_arma]];
            G_arima = [zeros(D,1); G_arma];
            scaling_param = 0.5;
            Q_arima = scaling_param * variance * (G_arima * G_arima.');
            H_arima = zeros(1,size(F_arima,1));
            H_arima(1:D+1) = 1;

            F = blkdiag(F_los,F_arima);
            Q = blkdiag(Q_los,Q_arima);
            H = [1, zeros(1,size(F_los,1) - 1), H_arima];
            R = diag(compute_phase_variances(general_config.C_over_N0_array_dBHz, sampling_interval));
            W = zeros(size(F,1),1);
        case 'rbf'
            % NOTE: For now, does nothing. This part is under development.
            % This case well never happen, considering that there is an
            % error in get_kalman_pll_config.m that is rasing to prevent it.
        case 'none'
            % For 'none', no augmentation is applied.
            % In this simple implementation, we leave the related matrices undefined.
            F = F_los;
            Q = Q_los;
            H = [1, zeros(1, size(F_los,1)-1)];
            R = diag(compute_phase_variances(general_config.C_over_N0_array_dBHz, sampling_interval));
            W = zeros(size(F_los,1), 1);
        otherwise
            error("MATLAB:UndefinedAugmentationModel","Invalid Augmentation model.")
    end

    % Store Results in Output Struct
    kalman_pll_config.(general_config.scintillation_training_data_config.scintillation_model) = struct( ...
        'F', F, ...
        'Q', Q, ...
        'H', H, ...
        'R', R, ...
        'W', W, ...
        'augmentation_model_initializer', general_config.augmentation_model_initializer ...
    );
end

function kalman_pll_config = build_kalman_pll_config(kalman_pll_config, ...
    scintillation_training_data_config, var_minimum_order, ...
    var_maximum_order, C_over_N0_array_dBHz, F_los, Q_los)
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
%   scintillation_training_data_config - Struct containing scintillation model settings.
%       For CSM, expected fields:
%           scintillation_model - Must be 'CSM'
%           S4                  - Scintillation index (0 <= S4 <= 1)
%           tau0                - Signal decorrelation time (positive scalar)
%           simulation_time     - Duration of simulation (positive scalar)
%           sampling_interval   - Sampling interval for the Kalman filter (positive scalar)
%
%       For TPPSM, expected fields:
%           scintillation_model - Must be 'TPPSM'
%           scenario            - A string specifying the scenario ('Weak', 'Moderate', 'Severe')
%           simulation_time     - Duration of simulation (positive scalar)
%           sampling_interval   - Sampling interval for the Kalman filter (positive scalar)
%           is_refractive_effects_removed - Boolean flag (true or false)
%
%   var_minimum_order     - Minimum VAR model order (scalar integer >= 1).
%   var_maximum_order     - Maximum VAR model order (scalar integer >= var_minimum_order).
%   C_over_N0_array_dBHz  - Array of carrier-to-noise density ratios in dB-Hz.
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

    % Validate inputs
    validateattributes(F_los, {'numeric'}, {'2d'}, mfilename, 'F_los');
    validateattributes(Q_los, {'numeric'}, {'2d'}, mfilename, 'Q_los');

    % Preprocess Training Data
    training_data = preprocess_training_data(scintillation_training_data_config);

    % Fit VAR Model (if applicable)
    if strcmp(scintillation_training_data_config.scintillation_model, 'none')
        intercept_vector = []; 
        var_coefficient_matrices = []; 
        var_covariance_matrices = [];
    else
        [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = ...
            arfit(training_data, var_minimum_order, var_maximum_order);
    end
    
    % Construct State Transition and Process Noise Matrices
    [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices( ...
        var_coefficient_matrices, var_covariance_matrices);

    % Construct Full Kalman Filter Matrices
    [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, ...
        intercept_vector, var_states_amount, var_model_order, ...
        C_over_N0_array_dBHz, scintillation_training_data_config.sampling_interval);

    % Store Results in Output Struct
    modelField = scintillation_training_data_config.scintillation_model;
    kalman_pll_config.(modelField) = struct( ...
        'F_los', F_los, ...
        'Q_los', Q_los, ...
        'F_var', F_var, ...
        'Q_var', Q_var, ...
        'F', F, ...
        'Q', Q, ...
        'H', H, ...
        'R', R, ...
        'W', W, ...
        'intercept_vector', intercept_vector, ...
        'var_model_order', var_model_order, ...
        'var_states_amount', var_states_amount ...
    );
end

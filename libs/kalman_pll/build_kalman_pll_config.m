function kalman_pll_config = build_kalman_pll_config(kalman_pll_config, ...
    scintillation_training_data_config, var_minimum_order, ...
    var_maximum_order, C_over_N0_array_dBHz, sampling_interval, ...
    F_los, Q_los)
    % build_kalman_pll_config
    % Computes the Kalman filter settings and VAR model matrices based on
    % the selected scintillation model and input configuration.
    %
    % Syntax:
    %   [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = ...
    %       build_kalman_pll_config(training_scint_model, scintillation_training_data_config, var_minimum_order, ...
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
    %   training_scint_model                      - Selected scintillation model ('CSM', 'TPPSM').
    %   scintillation_training_data_config - Cell array of scintillation model settings, such as:
    %                                         {S4, tau0, simulation_time, sampling_interval}.
    %   var_minimum_order                - Minimum VAR model order.
    %   var_maximum_order                - Maximum VAR model order.
    %   C_over_N0_array_dBHz             - Array of carrier-to-noise density ratios in dB-Hz.
    %   sampling_interval                - Sampling interval for the Kalman filter in seconds.
    %   F_los, Q_los                     - LOS dynamics state transition and covariance matrices.
    %   is_refractive_effects_removed    - Boolean flag indicating whether to remove refractive effects for TPPSM.
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
    %     and removes refractive effects if specified (for TPPSM only).
    %   - VAR model parameters are computed using the `arfit` function.
    %   - Full Kalman matrices are constructed by combining LOS and VAR model dynamics.
    %
    % Examples:
    %   % Compute settings for the 'CSM' model:
    %   [F_var, Q_var, F, Q, H, R, intercept_vector, var_states_amount, var_model_order] = ...
    %       build_kalman_pll_config('CSM', {0.8, 0.7, 300, 0.01}, 1, 6, [35], 0.01, F_los, Q_los, false);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    validateattributes(kalman_pll_config, {'struct'}, {'nonempty'}, mfilename, 'kalman_pll_config');
    
    validateattributes(scintillation_training_data_config, {'struct'}, {'nonempty'}, mfilename, 'scintillation_training_data_config');
    validateattributes(var_minimum_order, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'var_minimum_order');
    validateattributes(var_maximum_order, {'numeric'}, {'scalar', 'integer', '>=', var_minimum_order}, mfilename, 'var_maximum_order');
    validateattributes(C_over_N0_array_dBHz, {'numeric'}, {'vector', 'real', 'positive'}, mfilename, 'C_over_N0_array_dBHz');
    validateattributes(sampling_interval, {'numeric'}, {'scalar', 'real', 'positive'}, mfilename, 'sampling_interval');
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
        C_over_N0_array_dBHz, sampling_interval);

    % Store Results in Output Struct
    kalman_pll_config.(scintillation_training_data_config.scintillation_model) = struct( ...
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
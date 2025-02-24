function [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, var_states_amount, var_model_order, C_over_N0_array_dBHz, sampling_interval)
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
    W = [zeros(size(F_los, 1),1);intercept_vector;zeros(var_states_amount*(var_model_order-1),1)];
end
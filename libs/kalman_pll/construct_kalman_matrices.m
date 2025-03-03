function [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, var_states_amount, var_model_order, C_over_N0_array_dBHz, sampling_interval)
    % construct_kalman_matrices
    % Constructs the full Kalman filter matrices.
    %
    % Syntax:
    %   [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, var_states_amount, var_model_order, C_over_N0_array_dBHz, sampling_interval)
    %
    % Description:
    %   This function combines LOS dynamics and VAR model matrices to construct 
    %   the full state transition matrix (F), process noise covariance matrix (Q),
    %   measurement matrix (H), measurement noise covariance matrix (R), and an
    %   augmented intercept vector (W) for the Kalman filter.
    %
    % Inputs:
    %   F_los, Q_los     - LOS dynamics state transition and covariance matrices.
    %                      F_los must be a nonempty 2D numeric matrix.
    %                      Q_los must be a nonempty, square, symmetric, and positive semi-definite numeric matrix.
    %   F_var, Q_var     - VAR model state transition and covariance matrices.
    %                      F_var must be a nonempty 2D numeric matrix.
    %                      Q_var must be a nonempty, square, symmetric, and positive semi-definite numeric matrix.
    %   intercept_vector - Column vector of intercept terms.
    %   var_states_amount- Scalar integer (>=1) indicating the number of VAR states.
    %   var_model_order  - Scalar integer (>=1) indicating the VAR model order.
    %   C_over_N0_array_dBHz - Numeric vector (positive) of C/N0 values in dB-Hz.
    %   sampling_interval- Positive numeric scalar representing the sampling interval.
    %
    % Outputs:
    %   F - Full state transition matrix (blkdiag of F_los and F_var).
    %   Q - Full process noise covariance matrix (blkdiag of Q_los and Q_var).
    %   H - Measurement matrix constructed as:
    %         [1, zeros(1, size(F_los,1)-1), 1, zeros(1, var_states_amount*var_model_order-1)]
    %   R - Measurement noise covariance matrix computed as a diagonal matrix using
    %         compute_phase_variances(C_over_N0_array_dBHz, sampling_interval).
    %   W - Augmented intercept vector constructed as:
    %         [zeros(size(F_los,1),1); intercept_vector; zeros(var_states_amount*(var_model_order-1),1)]
    %
    % Example:
    %   [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_var, Q_var, intercept_vector, 3, 2, [35], 0.01);
    %
    % Author: Rodrigo de Lima Florindo
    % ORCID: https://orcid.org/0000-0003-0412-5583
    % Email: rdlfresearch@gmail.com

    % Validate inputs
    validateattributes(F_los, {'numeric'}, {'2d', 'nonempty'}, mfilename, 'F_los');
    validateattributes(Q_los, {'numeric'}, {'2d', 'nonempty', 'square'}, mfilename, 'Q_los');
    if ~isnan(var_model_order)
        validateattributes(F_var, {'numeric'}, {'2d', 'nonempty'}, mfilename, 'F_var');
        validateattributes(Q_var, {'numeric'}, {'2d', 'nonempty', 'square'}, mfilename, 'Q_var');
        validateattributes(intercept_vector, {'numeric'}, {'column'}, mfilename, 'intercept_vector');
        validateattributes(var_model_order, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'var_model_order');
        validateattributes(var_states_amount, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'var_states_amount');
    end
    
    % Ensure Q_los is symmetric and positive semi-definite
    if ~isequal(Q_los, Q_los') 
        error('construct_kalman_matrices:QlosNotSymmetric', 'Q_los must be symmetric.');
    end
    if any(eig(Q_los) < -1e-10) % Allow small numerical errors
        error('construct_kalman_matrices:QlosNotPositiveSemiDefinite', 'Q_los must be positive semi-definite.');
    end

    % Ensure Q_var is symmetric and positive semi-definite
    if ~isequal(Q_var, Q_var') 
        error('construct_kalman_matrices:QvarNotSymmetric', 'Q_var must be symmetric.');
    end
    if any(eig(Q_var) < -1e-10) 
        error('construct_kalman_matrices:QvarNotPositiveSemiDefinite', 'Q_var must be positive semi-definite.');
    end

    % Combine LOS and VAR state transitions
    F = blkdiag(F_los, F_var);
    
    % Combine LOS and VAR covariance matrices
    Q = blkdiag(Q_los, Q_var);
    
    % Compute measurement noise covariance matrix R from C/N0 values and sampling_interval.
    R = diag(compute_phase_variances(C_over_N0_array_dBHz, sampling_interval));

    if ~isnan(var_model_order)
        % Construct measurement matrix H.
        % H is defined as [1, zeros(1, size(F_los,1)-1), 1, zeros(1, var_states_amount*var_model_order-1)]
        H = [1, zeros(1, size(F_los,1)-1), 1, zeros(1, var_states_amount*var_model_order - 1)];
        
        % Construct the augmented intercept vector W.
        W = [zeros(size(F_los,1), 1); intercept_vector; zeros(var_states_amount*(var_model_order-1),1)];
    else
        H = [1, zeros(1, size(F_los,1)-1)];
        W = zeros(size(F_los,1), 1);
    end

    
end

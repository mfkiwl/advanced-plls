function [F_var, Q_var, var_states_amount, var_model_order] = build_var_matrices(var_coefficient_matrices, var_covariance_matrices)
    % build_var_matrices
    % Constructs VAR model matrices from coefficient and covariance matrices.
    %
    % Syntax:
    %   [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(var_coefficient_matrices, var_covariance_matrices)
    %
    % Description:
    %   This function generates the state transition matrix (F_var) and the process 
    %   noise covariance matrix (Q_var) for a VAR (Vector Autoregressive) model.
    %
    % Inputs:
    %   var_coefficient_matrices - Coefficient matrices for the VAR model. It should be a
    %                              nonempty 2D numeric matrix of size [n x (n * p)], where
    %                              n is the number of states and p is the VAR model order.
    %   var_covariance_matrices  - Covariance matrix for the VAR model. It should be a
    %                              nonempty square numeric matrix of size [n x n].
    %
    % Outputs:
    %   F_var           - Augmented state transition matrix for the VAR model.
    %   Q_var           - Augmented process noise covariance matrix for the VAR model.
    %   var_states_amount - Number of states (n) in the VAR model.
    %   var_model_order  - Order (p) of the VAR model.
    %
    % Notes:
    %   - The function constructs augmented matrices to handle higher-order VAR models.
    %   - Q_var includes zero padding for augmented states.
    %
    % Examples:
    %   [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(coeff_matrix, cov_matrix);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Validate inputs
    validateattributes(var_coefficient_matrices, {'numeric'}, {'2d'}, mfilename, 'var_coefficient_matrices');
    validateattributes(var_covariance_matrices, {'numeric'}, {'2d', 'square'}, mfilename, 'var_covariance_matrices');
    
    % Determine the number of states (n) from the coefficient matrix
    var_states_amount = size(var_coefficient_matrices, 1);
    
    % Check that the number of columns is an integer multiple of the number of states.
    numCols = size(var_coefficient_matrices, 2);
    if mod(numCols, var_states_amount) ~= 0
        error('build_var_matrices:InvalidDimensions', ...
            'The number of columns in var_coefficient_matrices must be an integer multiple of the number of rows.');
    end
    
    % Compute the VAR model order (p)
    var_model_order = numCols / var_states_amount;
    if ~isnan(var_model_order)
        validateattributes(var_model_order, {'numeric'}, {'scalar', 'integer', '>=', 1}, mfilename, 'var_model_order');
    end
    
    % Construct augmented state transition matrix
    % F_var = [A; [I_{n*(p-1)} , 0_{n*(p-1) x n}]]
    if var_model_order > 1
        F_var = [var_coefficient_matrices; [eye(var_states_amount * (var_model_order - 1)), zeros(var_states_amount * (var_model_order - 1), var_states_amount)]];
    else
        F_var = var_coefficient_matrices;  % For p == 1, no augmentation is needed.
    end

    % Construct augmented process noise covariance matrix
    % Q_var = blkdiag(Sigma, zeros(n*(p-1)))
    if var_model_order > 1
        Q_var = blkdiag(var_covariance_matrices, zeros(var_states_amount * (var_model_order - 1)));
    else
        Q_var = var_covariance_matrices;
    end
end

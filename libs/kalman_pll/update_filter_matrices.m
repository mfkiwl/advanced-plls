function [updated_F, updated_Q, updated_W] = update_filter_matrices(initial_F, initial_Q, initial_W, step, state_estimates, online_mdl_learning_cfg, augmentation_model_initializer)
% update_filter_matrices
%
% Syntax:
%   [updated_F, updated_Q, updated_W] = update_filter_matrices(initial_F, initial_Q, initial_W, step, state_estimates, online_mdl_learning_cfg, augmentation_model_initializer)
%
% Description:
%   Updates the state transition matrix (F), process noise covariance (Q), and intercept vector (W)
%   using online model learning. For AR models (identified by 'arfit' or 'aryule'), this function
%   supports two update methods based on the learning_method provided in online_mdl_learning_cfg:
%   'sliding_window' and 'block_window'. In the sliding_window method, AR coefficients and noise
%   variance are estimated using aryule. In the block_window method, arfit is used to estimate the
%   parameters and then construct VAR matrices. The update logic for RBF and Kalman models is under development.
%
% Inputs:
%   initial_F                - Initial state transition matrix.
%   initial_Q                - Initial state noise covariance matrix.
%   initial_W                - Initial intercept vector.
%   step                     - Current time step (integer).
%   state_estimates          - Matrix of state estimates (each row corresponds to a time step).
%   online_mdl_learning_cfg  - Struct with online learning configuration. For AR models, required fields are:
%                                learning_method - 'sliding_window' or 'block_window'.
%                                window_size     - Window size for the update.
%   augmentation_model_initializer - Struct with offline augmentation model configuration.
%                                Required fields:
%                                  id           - Augmentation model identifier ('arfit', 'aryule', 'rbf', 'kalman').
%                                  model_params - Struct containing model parameters (e.g., model_order).
%
% Outputs:
%   updated_F - Updated state transition matrix.
%   updated_Q - Updated process noise covariance matrix.
%   updated_W - Updated intercept vector.
%
% Notes:
%   - For AR models, the sliding_window method estimates AR coefficients via aryule, while the block_window method uses arfit
%     and a subsequent construction of VAR matrices via construct_var_matrices.
%   - The update logic for RBF and Kalman models is currently under development.
%
% Example:
%   % Example usage:
%   [upd_F, upd_Q, upd_W] = update_filter_matrices(F, Q, W, step, state_estimates, online_mdl_learning_cfg, augmentation_model_initializer);
%
% Dependencies:
%   aryule, arfit, construct_var_matrices, blkdiag.

    model_id = lower(char(augmentation_model_initializer.id));
    
    switch model_id
        case 'arfit'
            learning_method = lower(char(online_mdl_learning_cfg.learning_method));
            switch learning_method
                case 'sliding_window'
                if step >= online_mdl_learning_cfg.window_size
                    window_data = state_estimates(step - online_mdl_learning_cfg.window_size + 1 : step, :);
                    % Use arfit to estimate parameters for a block.
                    [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = ...
                        arfit(window_data, augmentation_model_initializer.model_params.model_order, augmentation_model_initializer.model_params.model_order);
                    
                    % Construct VAR state transition and process noise matrices.
                    [updated_F_var, updated_Q_var, var_states_amount] = construct_var_matrices(var_coefficient_matrices, var_covariance_matrices);
                    
                    wiener_model_order = size(initial_F, 2) - augmentation_model_initializer.model_params.model_order;
                    updated_F = blkdiag(initial_F(1:wiener_model_order, 1:wiener_model_order), updated_F_var);
                    updated_Q = blkdiag(initial_Q(1:wiener_model_order, 1:wiener_model_order), updated_Q_var);
                    
                    % Define var_model_order for clarity.
                    var_model_order = augmentation_model_initializer.model_params.model_order;
                    updated_W = [zeros(wiener_model_order, 1); intercept_vector; zeros(var_states_amount * (var_model_order - 1), 1)];
                else
                    updated_F = initial_F;
                    updated_Q = initial_Q;
                    updated_W = initial_W;
                end
                    
            case 'block_window'
                block_size = online_mdl_learning_cfg.window_size;
                if mod(step, block_size) == 0
                    window_data = state_estimates(step - block_size + 1 : step, :);
                    % Use arfit to estimate parameters for a block.
                    [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = ...
                        arfit(window_data, augmentation_model_initializer.model_params.model_order, augmentation_model_initializer.model_params.model_order);
                    
                    % Construct VAR state transition and process noise matrices.
                    [updated_F_var, updated_Q_var, var_states_amount] = construct_var_matrices(var_coefficient_matrices, var_covariance_matrices);
                    
                    wiener_model_order = size(initial_F, 2) - augmentation_model_initializer.model_params.model_order;
                    updated_F = blkdiag(initial_F(1:wiener_model_order, 1:wiener_model_order), updated_F_var);
                    updated_Q = blkdiag(initial_Q(1:wiener_model_order, 1:wiener_model_order), updated_Q_var);
                    
                    % Define var_model_order for clarity.
                    var_model_order = augmentation_model_initializer.model_params.model_order;
                    updated_W = [zeros(wiener_model_order, 1); intercept_vector; zeros(var_states_amount * (var_model_order - 1), 1)];
                else
                    updated_F = initial_F;
                    updated_Q = initial_Q;
                    updated_W = initial_W;
                end
            case 'kalman'
                % NOTE: Under development. Kalman-based update is not yet implemented.
                updated_F = initial_F;
                updated_Q = initial_Q;
                updated_W = initial_W;
            end
        case 'aryule'
            learning_method = lower(char(online_mdl_learning_cfg.learning_method));
            switch learning_method
                case 'sliding_window'
                    if step >= online_mdl_learning_cfg.window_size
                        % Determine the number of Wiener states.
                        wiener_model_order = size(initial_F, 2) - augmentation_model_initializer.model_params.model_order;
                        window_data = state_estimates(step - online_mdl_learning_cfg.window_size + 1 : step, wiener_model_order + 1);
                        % Estimate AR coefficients and noise variance using aryule.
                        [initial_coeffs, yule_variance] = aryule(window_data, augmentation_model_initializer.model_params.model_order);
                        % Format coefficients (exclude the first coefficient).
                        formatted_coeffs = -initial_coeffs(2:end);
                        
                        % Build AR state transition matrix.
                        updated_F_ar = [formatted_coeffs; [eye(augmentation_model_initializer.model_params.model_order - 1), zeros(augmentation_model_initializer.model_params.model_order - 1, 1)]];
                        % Build AR process noise covariance.
                        updated_Q_ar = zeros(augmentation_model_initializer.model_params.model_order);
                        updated_Q_ar(1,1) = yule_variance;
                        
                        % Combine the Wiener and AR parts.
                        updated_F = blkdiag(initial_F(1:wiener_model_order, 1:wiener_model_order), updated_F_ar);
                        updated_Q = blkdiag(initial_Q(1:wiener_model_order, 1:wiener_model_order), updated_Q_ar);
                        updated_W = initial_W;  % The intercept remains unchanged.
                        %step
                    else
                        updated_F = initial_F;
                        updated_Q = initial_Q;
                        updated_W = initial_W;
                    end
                    
                case 'block_window'
                    block_size = online_mdl_learning_cfg.window_size;
                    if mod(step, block_size) == 0
                        wiener_model_order = size(initial_F, 2) - augmentation_model_initializer.model_params.model_order;
                        window_data = state_estimates(step - block_size + 1 : step, wiener_model_order + 1);
                        % Estimate AR coefficients and noise variance using aryule.
                        [initial_coeffs, yule_variance] = aryule(window_data, augmentation_model_initializer.model_params.model_order);
                        % Format coefficients (exclude the first coefficient).
                        formatted_coeffs = -initial_coeffs(2:end);
                        
                        % Build AR state transition matrix.
                        updated_F_ar = [formatted_coeffs; [eye(augmentation_model_initializer.model_params.model_order - 1), zeros(augmentation_model_initializer.model_params.model_order - 1, 1)]];
                        % Build AR process noise covariance.
                        updated_Q_ar = zeros(augmentation_model_initializer.model_params.model_order);
                        updated_Q_ar(1,1) = yule_variance;
                        
                        updated_F = blkdiag(initial_F(1:wiener_model_order, 1:wiener_model_order), updated_F_ar);
                        updated_Q = blkdiag(initial_Q(1:wiener_model_order, 1:wiener_model_order), updated_Q_ar);
                        updated_W = initial_W;
                    else
                        updated_F = initial_F;
                        updated_Q = initial_Q;
                        updated_W = initial_W;
                    end
                otherwise
                    error('update_filter_matrices:invalid_learning_method', 'Invalid learning method for AR models: %s', learning_method);
            end
        case 'none'
            % Does nothing
            updated_F = initial_F;
            updated_Q = initial_Q;
            updated_W = initial_W;
        otherwise
            error('update_filter_matrices:invalid_model_id', 'Unsupported model id: %s', model_id);
    end
end
function aug_data = get_augmentation_model(training_data, general_config, sampling_interval)
% compute_augmentation computes augmentation matrices and related info based on the chosen model.
%
% Inputs:
%   training_data - Preprocessed training data.
%   general_config - Struct holding configuration parameters.
%   F_los, Q_los  - LOS dynamics matrices (computed from the discrete Wiener model).
%   sampling_interval - The sampling interval.
%
% Output:
%   aug_data - A struct containing augmentation matrices and additional parameters.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    switch lower(string(general_config.augmentation_model_initializer.id))
        case 'arfit'
            % Use ARFIT to estimate a VAR model.
            model_order = general_config.augmentation_model_initializer.model_params.model_order;
            [intercept_vector, var_coefficient_matrices, var_covariance_matrices] = ...
                arfit(training_data, model_order, model_order);
            [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(var_coefficient_matrices, var_covariance_matrices);
            aug_data.F_aug = F_var;
            aug_data.Q_aug = Q_var;
            aug_data.intercept = intercept_vector;
            % For ARFIT, we pass the computed intercept.
            aug_data.augInfo = struct('states_amount', var_states_amount, ...
                                     'model_order', var_model_order);
            aug_data.augmentation_type = 'arfit';

        case 'aryule'
            % Use ARYULE to estimate the VAR model parameters.
            model_order = general_config.augmentation_model_initializer.model_params.model_order;
            [ar_coefficients, ar_variance] = aryule(training_data, model_order);
            [F_var, Q_var, var_states_amount, var_model_order] = construct_var_matrices(...
                -ar_coefficients(2:end), ar_variance);
            aug_data.F_aug = F_var;
            aug_data.Q_aug = Q_var;
            aug_data.intercept = 0;
            % For ARYULE, set the intercept to zero.
            aug_data.augInfo = struct('states_amount', var_states_amount, ...
                                     'model_order', var_model_order);
            aug_data.augmentation_type = 'aryule';

        case 'kinematic'
            % Kinematic augmentation: use an alternative LOS-based Wiener model.
            L_aug = 1;  % single-frequency tracking
            M_aug = general_config.augmentation_model_initializer.model_params.wiener_mdl_order;
            sigma_aug = [0, 0, general_config.augmentation_model_initializer.model_params.process_noise_variance];
            delta_aug = 1;
            [F_wiener_aug, Q_wiener_aug] = get_discrete_wiener_model(L_aug, M_aug, sampling_interval, sigma_aug, delta_aug);
            aug_data.F_aug = F_wiener_aug;
            aug_data.Q_aug = Q_wiener_aug;
            aug_data.augmentation_type = 'kinematic';

        case 'arima'
            % ARIMA-based augmentation.
            p = general_config.augmentation_model_initializer.model_params.p;
            D = general_config.augmentation_model_initializer.model_params.D;
            q = general_config.augmentation_model_initializer.model_params.q;
            mdl = arima(p, D, q);
            mdl.Constant = 0;
            est_mdl = estimate(mdl, training_data);
            r = max(p, q + 1);
            ar_params = cell2mat(est_mdl.AR).';
            ma_params = cell2mat(est_mdl.MA).';
            variance = est_mdl.Variance;
            phi_vec = zeros(r,1);
            phi_vec(1:length(ar_params)) = ar_params;
            theta_vec = zeros(1, r-1);
            theta_vec(1:length(ma_params)) = ma_params;
            F_arma = [phi_vec, [eye(r-1); zeros(1, r-1)]];
            G_arma = [1, theta_vec].';
            % Build ARIMA state-space matrices:
            upper_right = zeros(D, r);
            upper_right(:,1) = 1;
            F_arima = [triu(ones(D)), upper_right; [zeros(r, D), F_arma]];
            G_arima = [zeros(D,1); G_arma];
            scaling_param = 0.5;
            Q_arima = scaling_param * variance * (G_arima * G_arima.');
            aug_data.F_aug = F_arima;
            aug_data.Q_aug = Q_arima;
            % Store the extra measurement vector used in the ARIMA branch.
            H_arima = zeros(1, size(F_arima,1));
            H_arima(1:D+1) = 1;
            aug_data.H_aug = H_arima;
            aug_data.augmentation_type = 'arima';

        case 'rbf'
            % RBF branch is not yet implemented.
            error('MATLAB:RBFUnavailable', 'RBF model initializer is still under development.');

        case 'none'
            % No augmentation: return an empty augmentation structure (or flag it as "none").
            aug_data = struct();
            aug_data.augmentation_type = 'none';

        otherwise
            error('MATLAB:UndefinedAugmentationModel', 'Invalid augmentation model: %s', general_config.augmentation_model_initializer.id);
    end
end

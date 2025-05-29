function initial_estimates = get_initial_estimates(general_config, kalman_pll_config)
% get_initial_estimates
% Generates initial state estimates and covariance matrices for Kalman filter-based PLL.
%
% Syntax:
%   initial_estimates = get_initial_estimates(general_config, kalman_pll_config)
%
% Description:
%   This function generates initial estimates for the state vector and covariance matrix
%   used by the Kalman filter. It supports two modes:
%     1. Generating random initial estimates based on uniform distributions.
%     2. Using perfect estimates for the initial state.
%   LOS dynamics variances are computed from the uniform distributions specified in
%   initial_states_distributions_boundaries. The variance for VAR phase states is assumed
%   to be (pi^2/3), corresponding to a uniform distribution in [-pi, pi].
%
% Inputs:
%   general_config - Struct containing configuration details with the following fields:
%       - is_generate_random_initial_estimates: Boolean flag indicating whether to generate random
%         initial estimates based on initial_states_distributions_boundaries.
%       - initial_states_distributions_boundaries: Cell array of boundaries for uniform distributions 
%         to generate random initial estimates of the LOS states. Each element must be a 1×2 numeric
%         vector with the first element less than the second.
%       - expected_doppler_profile: Non-empty numeric vector.
%       - scintillation_training_data_config: Struct with field 'scintillation_model'
%         (e.g., 'CSM', 'TPPSM', or 'none').
%
%   kalman_pll_config - Struct containing Kalman filter settings, which must include:
%       - A field corresponding to the training scintillation model (e.g., 'CSM')
%         that contains the following fields:
%           var_states_amount: Number of VAR model states.
%           var_model_order: Order of the VAR model.
%
% Outputs:
%   initial_estimates - Struct containing:
%       * x_hat_init: Initial state vector estimate (column vector). Its length equals:
%                     length(expected_doppler_profile) + (var_states_amount * var_model_order).
%       * P_hat_init: Initial covariance matrix of size [L+N x L+N], where L is the length of
%                     expected_doppler_profile and N = var_states_amount * var_model_order.
%
% References:
%  [1] Durbin, James & Koopman, Siem Jan. (2001). Time Series Analysis 
%      by State Space Methods. 
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Validate each boundary.
    for i = 1:length(general_config.initial_states_distributions_boundaries)
        b = general_config.initial_states_distributions_boundaries{i};
        validateattributes(b, {'numeric'}, {'vector', 'numel', 2}, mfilename, 'initial_states_distributions_boundaries');
        if b(1) >= b(2)
            error('get_initial_estimates:InvalidBoundaries', ...
                'Each boundary must have its first element less than its second.');
        end
    end

    % Validate kalman_pll_config is nonempty.
    validateattributes(kalman_pll_config, {'struct'}, {'nonempty'}, mfilename, 'kalman_pll_config');
    
    % Ensure that training_scint_model is a field in kalman_pll_config.
    model_name = char(general_config.scintillation_training_data_config.scintillation_model);
    if ~isfield(kalman_pll_config, model_name)
        error('get_initial_estimates:MissingScintModelField', ...
            'kalman_pll_config does not contain the model field "%s".', model_name);
    end
    
    % Initialize output struct.
    initial_estimates = struct('x_hat_init', [], 'P_hat_init', []);

    L = numel(general_config.expected_doppler_profile);
    switch general_config.augmentation_model_initializer.id
        case {'arfit', 'aryule'}
            % Determine lengths.
            var_states = general_config.discrete_wiener_model_config{1};
            var_order = general_config.augmentation_model_initializer.model_params.model_order;
            N = var_states * var_order;
            if general_config.is_generate_random_initial_estimates
                % Generate random errors for LOS states based on uniform distributions.
                random_errors = zeros(L,1);
                for i = 1:L
                    bound = general_config.initial_states_distributions_boundaries{i};
                    random_errors(i) = unifrnd(bound(1), bound(2));
                end
                
                % x_hat_init: LOS states adjusted by random error, and VAR states set to zero.
                initial_estimates.x_hat_init = [general_config.expected_doppler_profile.' - random_errors; zeros(N,1)];
                
                % LOS variances computed from uniform distribution variance: (b2 - b1)^2 / 12.
                los_variances = cellfun(@(b) (b(2) - b(1))^2 / 12, general_config.initial_states_distributions_boundaries);
                % For VAR states, variance is assumed to be (pi^2/3).
                kinematic_variance = (pi^2/3);
                initial_estimates.P_hat_init = blkdiag(diag(los_variances), kinematic_variance * eye(N));
            else
                % Use perfect estimates: LOS states are exactly the real doppler profile.
                initial_estimates.x_hat_init = [general_config.expected_doppler_profile.'; zeros(N,1)];
                
                % Covariance: Perfect knowledge for LOS states (except first state has variance pi^2/3),
                % VAR states get variance pi^2/3.
                LOS_cov = diag([pi^2/3, zeros(1, L-1)]);
                initial_estimates.P_hat_init = blkdiag(LOS_cov, (pi^2/3) * eye(N));
            end

        case 'rbf'
            % NOTE: This part is under development.
        case 'kinematic'
            if general_config.is_generate_random_initial_estimates
                % Generate random errors for LOS states based on uniform distributions.
                random_errors = zeros(L,1);
                for i = 1:L
                    bound = general_config.initial_states_distributions_boundaries{i};
                    random_errors(i) = unifrnd(bound(1), bound(2));
                end
                
                % x_hat_init: LOS states adjusted by random error, and VAR states set to zero.
                initial_estimates.x_hat_init = [general_config.expected_doppler_profile.' - random_errors; zeros(general_config.augmentation_model_initializer.model_params.wiener_mdl_order,1)];
                
                % LOS variances computed from uniform distribution variance: (b2 - b1)^2 / 12.
                los_variances = cellfun(@(b) (b(2) - b(1))^2 / 12, general_config.initial_states_distributions_boundaries);
                % For VAR states, variance is assumed to be (pi^2/3).
                initial_estimates.P_hat_init = blkdiag(diag(los_variances), zeros(general_config.augmentation_model_initializer.model_params.wiener_mdl_order));
            else
                % Use perfect estimates: LOS states are exactly the real doppler profile.
                initial_estimates.x_hat_init = [general_config.expected_doppler_profile.'; zeros(general_config.augmentation_model_initializer.model_params.wiener_mdl_order,1)];
                
                % Covariance: Perfect knowledge for LOS states (except first state has variance pi^2/3),
                % VAR states get variance pi^2/3.
                LOS_cov = diag([pi^2/3, zeros(1, L-1)]);
                initial_estimates.P_hat_init = blkdiag(LOS_cov, zeros(general_config.augmentation_model_initializer.model_params.wiener_mdl_order));
            end
        case 'arima'
            % ARIMA model orders that affect the dimensions of the
            % state vector: p and D. The total amount of lines is given
            % by p + D.

            % The development of the initial estimates and covariance
            % matrices was inspired by what is proposed in [Section 5.6.3, 1]
            p = general_config.augmentation_model_initializer.model_params.p;
            D = general_config.augmentation_model_initializer.model_params.D;
            if general_config.is_generate_random_initial_estimates
                % Generate random errors for LOS states based on uniform distributions.
                random_errors = zeros(L,1);
                for i = 1:L
                    bound = general_config.initial_states_distributions_boundaries{i};
                    random_errors(i) = unifrnd(bound(1), bound(2));
                end
                
                % x_hat_init: LOS states adjusted by random error, and VAR states set to zero.
                initial_estimates.x_hat_init = [general_config.expected_doppler_profile.' - random_errors; zeros(p + D,1)];
                
                % LOS variances computed from uniform distribution variance: (b2 - b1)^2 / 12.
                los_variances = cellfun(@(b) (b(2) - b(1))^2 / 12, general_config.initial_states_distributions_boundaries);
                
                kappa = 1e-2;
                P_inf = zeros(p + D);
                P_inf(1,1) = 1;

                Q = kalman_pll_config.(general_config.scintillation_training_data_config.scintillation_model).Q;
                Q_arima = Q(L+1 : end, L+1 : end);
                P_asterisk = Q_arima;

                P_arima_init = kappa*P_inf + P_asterisk;

                initial_estimates.P_hat_init = blkdiag(diag(los_variances), P_arima_init);
            else
                % Use perfect estimates: LOS states are exactly the real doppler profile.
                initial_estimates.x_hat_init = [general_config.expected_doppler_profile.'; zeros(p + D,1)];
                
                kappa = 1e-2;
                P_inf = zeros(p+D);
                P_inf(1,1) = 1;

                Q = kalman_pll_config.(general_config.scintillation_training_data_config.scintillation_model).Q;
                Q_arima = Q(L+1:end, L+1:end);
                P_asterisk = Q_arima;

                P_arima_init = kappa*P_inf + P_asterisk;

                % Covariance: Perfect knowledge for LOS states (except first state has variance pi^2/3),
                % VAR states get variance pi^2/3.
                LOS_cov = diag([pi^2/3, zeros(1, L-1)]);
                initial_estimates.P_hat_init = blkdiag(LOS_cov, P_arima_init);
            end
        case 'none'
            if general_config.is_generate_random_initial_estimates
                % Generate random errors for LOS states based on uniform distributions.
                random_errors = zeros(L,1);
                for i = 1:L
                    bound = general_config.initial_states_distributions_boundaries{i};
                    random_errors(i) = unifrnd(bound(1), bound(2));
                end
                
                % x_hat_init: LOS states adjusted by random error, and VAR states set to zero.
                initial_estimates.x_hat_init = general_config.expected_doppler_profile.' - random_errors;
                
                % LOS variances computed from uniform distribution variance: (b2 - b1)^2 / 12.
                los_variances = cellfun(@(b) (b(2) - b(1))^2 / 12, general_config.initial_states_distributions_boundaries);
                initial_estimates.P_hat_init = diag(los_variances);
            else
                % Use perfect estimates: LOS states are exactly the real doppler profile.
                initial_estimates.x_hat_init = general_config.expected_doppler_profile.';
                
                % Covariance: Perfect knowledge for LOS states (except first state has variance pi^2/3),
                LOS_cov = diag([pi^2/3, zeros(1, L-1)]);
                initial_estimates.P_hat_init = LOS_cov;
            end
            
    end
end

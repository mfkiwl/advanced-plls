function initial_estimates = handle_initial_estimates(config, kalman_pll_config)
    % handle_initial_estimates
    % Generates initial state estimates and covariance matrices for Kalman filter-based PLL.
    %
    % Syntax:
    %   initial_estimates = handle_initial_estimates(config, kalman_pll_config)
    %
    % Description:
    %   This function generates initial estimates for the state vector and covariance matrix
    %   used by the Kalman filter. It supports two modes: 
    %   1. Generating random initial estimates based on uniform distributions.
    %   2. Using perfect estimates for the initial state.
    %   The function accounts for line-of-sight (LOS) dynamics and VAR (Vector Autoregressive) states,
    %   with appropriate variances based on the configuration.
    %
    % Inputs:
    %   config - Struct containing configuration details with the following fields:
    %       - is_generate_random_initial_estimates: Boolean flag indicating whether to generate random
    %         initial estimates based on `initial_states_distributions_boundaries`.
    %       - initial_states_distributions_boundaries: Cell array of boundaries for uniform distributions 
    %         to generate random initial estimates of the Doppler profile. Example:
    %           {{[-pi, pi]}, {[-5, 5]}, {[-0.1, 0.1]}, {[-0.001, 0.001]}}
    %       - real_doppler_profile: Vector containing the real Doppler profile values.
    %       - training_scint_model: String specifying the scintillation model ('CSM', 'TPPSM', or 'none').
    %
    %   kalman_pll_config - Struct containing Kalman filter settings, which must include:
    %       - var_states_amount: Number of VAR model states.
    %       - var_model_order: Order of the VAR model.
    %
    % Outputs:
    %   initial_estimates - Struct containing:
    %       * x_hat_init: Initial state vector estimate, with LOS dynamics states followed by VAR states.
    %       * P_hat_init: Initial covariance matrix, with LOS variances and uniform variance for VAR states.
    %
    % Notes:
    %   - LOS dynamics variances are computed using the uniform distributions specified in
    %     `initial_states_distributions_boundaries`.
    %   - The variance for VAR phase states is assumed to be \((\pi^2 / 3)\), corresponding to a uniform
    %     distribution within \([- \pi, \pi]\).
    %
    % Examples:
    %   % Example configuration for initial estimates:
    %   config = struct( ...
    %       'is_generate_random_initial_estimates', true, ...
    %       'initial_states_distributions_boundaries', {{[-pi, pi], [-5, 5], [-0.1, 0.1]}}, ...
    %       'real_doppler_profile', [0.01, -0.02, 0.03], ...
    %       'training_scint_model', 'CSM' ...
    %   );
    %   kalman_pll_config = struct( ...
    %       'CSM', struct('var_states_amount', 2, 'var_model_order', 1) ...
    %   );
    %   initial_estimates = handle_initial_estimates(config, kalman_pll_config);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com
    
    initial_estimates = struct('x_hat_init', [], 'P_hat_init', []);
    if config.is_generate_random_initial_estimates
        %%%% Generate an error vector for the initial state estimates.
        random_doppler_profile_error = arrayfun(@(col) ...
            unifrnd(config.initial_states_distributions_boundaries{1, col}(1), ...
                    config.initial_states_distributions_boundaries{1, col}(2)), ...
            1:length(config.real_doppler_profile));

        %%% Create an initial estimates vector where the first elements
        % corresponds to the line-of-sight dynamics states and the other
        % ones to the VAR states. The VAR initial states are assumed to be
        % zero arbitrarily. Nevertheless, the initial values for the VAR
        % phase estimates are not very much important, since they can be
        % any number within the range of [-pi,pi]. If we want to refactor
        % this code to be used by an extended Kalman filter instead, we
        % need to refactor this part to comprehend the amplitude estimates.
        initial_estimates.x_hat_init = [ ...
            config.real_doppler_profile.' - random_doppler_profile_error.'; ...
            zeros(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order, 1) ...
        ];

        %%% Compute line-of-sight dynamics variances according to the
        % distributions from
        % `config.initial_states_distributions_boundaries`.
        los_variances = ...
        cellfun( ...
            @(nested) ...
            (nested(2) - nested(1))^2 / 12, ...
            config.initial_states_distributions_boundaries ...
        );
        
        %%% Construct diagonal covariance matrix
        % All VAR phase states are assumed to present a variance of (pi^2/3),
        % which refers to the variance of a uniform distribution that
        % lies within the range [-pi,pi]. 
        initial_estimates.P_hat_init = ...
        blkdiag( ...
            diag(los_variances), ...
            (pi^2/3) * eye(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order) ...
        );
    else
        % Build x_hat_init with the perfect estimates for the
        % line-of-sight dynamics
        initial_estimates.x_hat_init = [ ...
            config.real_doppler_profile.'; ...
            zeros(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order, 1)
        ];
        % Build P_hat_init considering that the doppler profile is
        % perfectly known, i.e., they have covariances elements equal to 
        % zero, with the exception of the actual line-of-sight and VAR phases.
        initial_estimates.P_hat_init = blkdiag( ...
            diag([(pi^2/3),zeros(1,length(config.real_doppler_profile) - 1)]), ...
            (pi^2/3) * eye(kalman_pll_config.(config.training_scint_model).var_states_amount * ...
            kalman_pll_config.(config.training_scint_model).var_model_order) ...
        );
    end
end
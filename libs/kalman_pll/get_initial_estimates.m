function initial_estimates = get_initial_estimates(config, kalman_pll_config)
% get_initial_estimates
% Generates initial state estimates and covariance matrices for Kalman filter-based PLL.
%
% Syntax:
%   initial_estimates = get_initial_estimates(config, kalman_pll_config)
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
%   initial_estimates = get_initial_estimates(config, kalman_pll_config);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 ORCID: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

    % Validate config is a nonempty struct
    validateattributes(config, {'struct'}, {'nonempty'}, mfilename, 'config');
    
    % Check required fields in config
    required_config_fields = {
        'is_generate_random_initial_estimates', ...
        'initial_states_distributions_boundaries', ...
        'real_doppler_profile', ...
        'training_scint_model'
    };
    for iField = 1:numel(required_config_fields)
        if ~isfield(config, required_config_fields{iField})
            error('get_initial_estimates:MissingConfigField', ...
                'The config struct must contain the field "%s".', required_config_fields{iField});
        end
    end

    % Validate types and sizes of config fields
    validateattributes(config.is_generate_random_initial_estimates, ...
        {'logical'}, {'scalar'}, mfilename, 'is_generate_random_initial_estimates');
    validateattributes(config.initial_states_distributions_boundaries, ...
        {'cell'}, {'nonempty'}, mfilename, 'initial_states_distributions_boundaries');
    validateattributes(config.real_doppler_profile, ...
        {'numeric'}, {'vector', 'finite', 'nonnan'}, mfilename, 'real_doppler_profile');
    validateattributes(config.training_scint_model, ...
        {'char','string'}, {'nonempty'}, mfilename, 'training_scint_model');

    % Check the number of boundaries matches the length of real_doppler_profile
    if numel(config.initial_states_distributions_boundaries) ~= numel(config.real_doppler_profile)
        error('get_initial_estimates:BoundaryProfileMismatch', ...
            'Number of boundaries must match length of real_doppler_profile.');
    end
    
    % Validate kalman_pll_config is nonempty
    validateattributes(kalman_pll_config, {'struct'}, {'nonempty'}, mfilename, 'kalman_pll_config');
    % Ensure training_scint_model is a field in kalman_pll_config
    model_name = char(config.training_scint_model);
    if ~isfield(kalman_pll_config, model_name)
        error('get_initial_estimates:MissingScintModelField', ...
            'kalman_pll_config does not contain the model field "%s".', model_name);
    end
    
    % Check that var_states_amount and var_model_order exist
    required_kalman_fields = {'var_states_amount','var_model_order'};
    for iField = 1:numel(required_kalman_fields)
        if ~isfield(kalman_pll_config.(model_name), required_kalman_fields{iField})
            error('get_initial_estimates:MissingKalmanField', ...
                'kalman_pll_config.%s is missing the field "%s".', ...
                model_name, required_kalman_fields{iField});
        end
    end

    initial_estimates = struct('x_hat_init', [], 'P_hat_init', []);

    if config.is_generate_random_initial_estimates
        % Generate error vector for the initial LOS states
        random_doppler_profile_error = arrayfun(@(col) ...
            unifrnd(config.initial_states_distributions_boundaries{1, col}(1), ...
                    config.initial_states_distributions_boundaries{1, col}(2)), ...
            1:length(config.real_doppler_profile));
        
        % Combine LOS states with zeroed VAR states in x_hat_init
        initial_estimates.x_hat_init = [ ...
            config.real_doppler_profile.' - random_doppler_profile_error.'; ...
            zeros(kalman_pll_config.(model_name).var_states_amount * ...
                   kalman_pll_config.(model_name).var_model_order, 1) ...
        ];
        
        % Compute LOS variances from the uniform distributions
        los_variances = cellfun(@(bound) (bound(2) - bound(1))^2 / 12, ...
            config.initial_states_distributions_boundaries);

        % Build covariance matrix: LOS dynamics + uniform for VAR states
        initial_estimates.P_hat_init = blkdiag( ...
            diag(los_variances), ...
            (pi^2/3) * eye(kalman_pll_config.(model_name).var_states_amount * ...
            kalman_pll_config.(model_name).var_model_order) ...
        );
    else
        % Use perfect estimates for LOS states
        initial_estimates.x_hat_init = [ ...
            config.real_doppler_profile.'; ...
            zeros(kalman_pll_config.(model_name).var_states_amount * ...
                   kalman_pll_config.(model_name).var_model_order, 1) ...
        ];

        % Build covariance matrix with zero for LOS states (perfect knowledge),
        % except for the first LOS state which we assume pi^2/3 variance, 
        % plus uniform for VAR states.
        initial_estimates.P_hat_init = blkdiag( ...
            diag([(pi^2/3), zeros(1, length(config.real_doppler_profile) - 1)]), ...
            (pi^2/3) * eye(kalman_pll_config.(model_name).var_states_amount * ...
            kalman_pll_config.(model_name).var_model_order) ...
        );
    end
end

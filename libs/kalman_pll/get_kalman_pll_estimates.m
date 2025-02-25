function [state_estimates, error_covariance_estimates] = get_kalman_pll_estimates(received_signal, kalman_pll_config, initial_estimates, training_scint_model)
    % get_kalman_pll_estimates
    % Generates Kalman filter state and error covariance estimates based on the
    % received signal and the provided Kalman PLL configuration.
    %
    % Syntax:
    %   [state_estimates, error_covariance_estimates] = get_kalman_pll_estimates(received_signal, kalman_pll_config, initial_estimates, training_scint_model)
    %
    % Description:
    %   This function applies the Kalman filter algorithm iteratively to the
    %   received signal to produce state estimates and error covariance estimates.
    %   It uses the configuration settings (F, Q, H, R, W) stored in the field of
    %   kalman_pll_config corresponding to training_scint_model.
    %
    % Inputs:
    %   received_signal     - Numeric 2D array (NxM) representing the received signal.
    %                         The function uses only the first column.
    %   kalman_pll_config   - Struct containing the Kalman filter configuration.
    %                         Must have a field named after training_scint_model.
    %   initial_estimates   - Struct with fields:
    %                             x_hat_init: Numeric column vector (initial state estimate)
    %                             P_hat_init: Numeric square matrix (initial error covariance)
    %
    % Outputs:
    %   state_estimates           - Numeric matrix of state estimates. Each row corresponds
    %                               to a time step; number of columns equals the length of x_hat_init.
    %   error_covariance_estimates- 3D numeric array containing error covariance estimates.
    %                               Its first dimension is time, and its second and third dimensions
    %                               match the size of P_hat_init.
    %
    % References: 
    % [1] Vilà-Valls J, Closas P, Curran JT. 2017. Multi-frequency GNSS robust carrier tracking 
    %     for ionospheric scintillation mitigation. J. Space Weather Space Clim. 7: A26
    % [2] Florindo, Rodrigo de Lima, Antreich, Felix, "Multi-Frequency Kalman Filter Carrier 
    %     Phase Tracking for Ionospheric Scintillation Mitigation and Monitoring," Proceedings
    %     of the 37th International Technical Meeting of the Satellite Division of The Institute
    %     of Navigation (ION GNSS+ 2024), Baltimore, Maryland, September 2024, pp. 3611-3625.
    %     https://doi.org/10.33012/2024.19899
    % Author: Rodrigo de Lima Florindo
    % ORCID: https://orcid.org/0000-0003-0412-5583
    % Email: rdlfresearch@gmail.com

    % Validate inputs
    validateattributes(received_signal, {'numeric'}, {'nonempty','2d'}, mfilename, 'received_signal');
    validateattributes(kalman_pll_config, {'struct'}, {'nonempty'}, mfilename, 'kalman_pll_config');
    validateattributes(initial_estimates, {'struct'}, {'nonempty'}, mfilename, 'initial_estimates');
    
    if ~isfield(initial_estimates, 'x_hat_init') || ~isfield(initial_estimates, 'P_hat_init')
        error('get_kalman_pll_estimates:MissingField', 'initial_estimates must have fields x_hat_init and P_hat_init.');
    end
    validateattributes(initial_estimates.x_hat_init, {'numeric'}, {'nonempty','vector'}, mfilename, 'initial_estimates.x_hat_init');
    % Ensure P_hat_init is a square matrix whose size matches the length of x_hat_init.
    n = numel(initial_estimates.x_hat_init);
    validateattributes(initial_estimates.P_hat_init, {'numeric'}, {'nonempty','2d','size',[n n]}, mfilename, 'initial_estimates.P_hat_init');
    
    if ~(ischar(training_scint_model) || isstring(training_scint_model))
        error('get_kalman_pll_estimates:InvalidType', 'training_scint_model must be a char or string.');
    end
    training_scint_model = char(training_scint_model); % Convert to char if needed.
    if ~isfield(kalman_pll_config, training_scint_model)
        error('get_kalman_pll_estimates:InvalidTrainingModel', 'kalman_pll_config does not have the field for training_scint_model: %s', training_scint_model);
    end

    % Retrieve Kalman filter matrices from configuration.
    configStruct = kalman_pll_config.(training_scint_model);
    requiredConfigFields = {'F','Q','H','R','W'};
    for i = 1:length(requiredConfigFields)
        if ~isfield(configStruct, requiredConfigFields{i})
            error('get_kalman_pll_estimates:MissingConfigField', ...
                'kalman_pll_config.%s is missing the field %s.', training_scint_model, requiredConfigFields{i});
        end
    end
    F = configStruct.F;
    Q = configStruct.Q;
    H = configStruct.H;
    R = configStruct.R;
    W = configStruct.W;
    
    % Preallocate output arrays.
    N = size(received_signal,1);
    state_estimates = zeros(N, numel(initial_estimates.x_hat_init));
    error_covariance_estimates = zeros(N, size(initial_estimates.P_hat_init,1), size(initial_estimates.P_hat_init,2));
    
    % Initialize estimates.
    x_hat_project_ahead = initial_estimates.x_hat_init;
    P_hat_project_ahead = initial_estimates.P_hat_init;
    
    % Loop through each time step.
    for step = 1:N
        if step > 1
            % Compute Kalman gain.
            K = P_hat_project_ahead * H.' * ((H * P_hat_project_ahead * H.' + R) \ eye(size(R,1)));
            % Update state estimate using the phase error from the received signal.
            x_hat_update = x_hat_project_ahead + K * angle(received_signal(step-1,1) * exp(-1j * (H * x_hat_project_ahead)));
            % Update error covariance.
            P_hat_update = P_hat_project_ahead - K * H * P_hat_project_ahead;
        else
            x_hat_update = x_hat_project_ahead;
            P_hat_update = P_hat_project_ahead;
        end
        % Project ahead.
        % NOTE: W is being used here to introduce the bias of the VAR
        % model. See [1,Section 3.3] and [2, Equation 20].
        x_hat_project_ahead = F * x_hat_update + W;
        P_hat_project_ahead = F * P_hat_update * F.' + Q;
        % Save estimates.
        state_estimates(step,:) = x_hat_project_ahead.';
        error_covariance_estimates(step,:,:) = P_hat_project_ahead;
    end
end

function [state_estimates, error_covariance_estimates] = get_kalman_pll_estimates(received_signal, kalman_pll_config, initial_estimates, training_scint_model, adaptive_config)
% get_kalman_pll_estimates
% Generates Kalman filter state and error covariance estimates based on the
% received signal and the provided Kalman PLL configuration.
%
% Syntax:
%   [state_estimates, error_covariance_estimates] = get_kalman_pll_estimates(received_signal, kalman_pll_config, initial_estimates, training_scint_model, adaptive_config)
%
% Description:
%   This function applies the Kalman filter algorithm iteratively to the
%   received signal to produce state estimates and error covariance estimates.
%   It uses the configuration settings (F, Q, H, R, W) stored in the field of
%   kalman_pll_config corresponding to the given training_scint_model.
%
%   The training_scint_model parameter can be one of:
%     'CSM'      - Carrier Scintillation Model with AR augmentation.
%     'TPPSM'    - TPPSM model with AR augmentation.
%     'none'     - Standard KF (no AR model augmentation).
%
%   In addition, an adaptive configuration is applied based on the provided
%   adaptive_config struct. The allowed options for adaptive_config.algorithm are:
%      'simplified' - Applies a simplified adaptive update (see placeholder below).
%      'NWBP'       - Not implemented (an error is thrown if chosen).
%      'none'       - No adaptive update is applied (i.e. R is used directly).
%
%   Note: If the 'simplified' option is chosen, and adaptive_config.algorithm is
%         not 'none', then the fields L1_C_over_N0_dBHz, sampling_interval, and threshold
%         must be provided.
%
% Inputs:
%   received_signal      - Numeric 2D array (NxM) representing the received signal.
%                          (Only the first column is used.)
%   kalman_pll_config    - Struct containing the Kalman filter configuration.
%                          Must have a field named after training_scint_model.
%   initial_estimates    - Struct with fields:
%                              x_hat_init: Numeric column vector (initial state estimate)
%                              P_hat_init: Numeric square matrix (initial error covariance)
%   training_scint_model - A string with the name of the scintillation model used to
%                          train the AR model. Options: {'CSM', 'TPPSM', 'none'}.
%   adaptive_config      - Struct with adaptive configuration options for the
%                          Kalman filter algorithm. It must have the fields:
%         .algorithm    - {'simplified','NWBP','none'}. If 'simplified', a simplified adaptive
%                         update is applied (using a placeholder based on [3, Equation 64]).
%                         If 'NWBP' is selected, an error is thrown (not yet implemented).
%                         If 'none' is selected, no adaptation is performed.
%         .hard_limited - Logical; if true, a hard-limited constraint (as in [5, Eq.24])
%                         is applied (using a placeholder large value in case of too-small noise).
%                         When algorithm is 'simplified', adaptive_config must also have:
%                            - L1_C_over_N0_dBHz
%                            - sampling_interval
%                            - threshold
%
% Outputs:
%   state_estimates            - Numeric matrix of state estimates. Each row corresponds
%                                to a time step; the number of columns equals the length of x_hat_init.
%   error_covariance_estimates - 3D numeric array containing error covariance estimates.
%                                Its first dimension is time, and the second and third dimensions
%                                match the size of P_hat_init.
%
% References: 
% [1] Vilà-Valls J, Closas P, Curran JT. 2017. Multi-frequency GNSS robust carrier tracking 
%     for ionospheric scintillation mitigation. J. Space Weather Space Clim. 7: A26
% [2] Florindo, Rodrigo de Lima, Antreich, Felix, "Multi-Frequency Kalman Filter Carrier 
%     Phase Tracking for Ionospheric Scintillation Mitigation and Monitoring," Proceedings
%     of the 37th International Technical Meeting of the Satellite Division of The Institute
%     of Navigation (ION GNSS+ 2024), Baltimore, Maryland, September 2024, pp. 3611-3625.
%     https://doi.org/10.33012/2024.19899
% [3] R. A. M. Lopes, F. Antreich, F. Fohlmeister, M. Kriegel and H. K. Kuga, "Ionospheric 
%     Scintillation Mitigation With Kalman PLLs Employing Radial Basis Function Networks," 
%     in IEEE Transactions on Aerospace and Electronic Systems, vol. 59, no. 5, pp. 6878-6893,
%     Oct. 2023, doi: 10.1109/TAES.2023.3281431
% [4] Locubirche-Serra, Sergi. Robust Carrier Tracking Techniques for GNSS Receivers affected by
%     Ionospheric Scintillation. Doctoral Thesis. https://ddd.uab.cat/record/235075?ln=en
% [5] Locubiche-Serra, S., Seco-Granados, G. & López-Salcedo, J.A. Performance assessment of a 
%     low-complexity autoregressive Kalman filter for GNSS carrier tracking using real 
%     scintillation time series. GPS Solut 26, 17 (2022). https://doi.org/10.1007/s10291-021-01193-0
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    %% Validate inputs
    validateattributes(received_signal, {'numeric'}, {'nonempty','2d'}, mfilename, 'received_signal');
    validateattributes(kalman_pll_config, {'struct'}, {'nonempty'}, mfilename, 'kalman_pll_config');
    validateattributes(initial_estimates, {'struct'}, {'nonempty'}, mfilename, 'initial_estimates');
    
    if ~isfield(initial_estimates, 'x_hat_init') || ~isfield(initial_estimates, 'P_hat_init')
        error('get_kalman_pll_estimates:MissingField', 'initial_estimates must have fields x_hat_init and P_hat_init.');
    end
    validateattributes(initial_estimates.x_hat_init, {'numeric'}, {'nonempty','vector'}, mfilename, 'initial_estimates.x_hat_init');
    n = numel(initial_estimates.x_hat_init);
    validateattributes(initial_estimates.P_hat_init, {'numeric'}, {'nonempty','2d','size',[n n]}, mfilename, 'initial_estimates.P_hat_init');
    
    if ~(ischar(training_scint_model) || isstring(training_scint_model))
        error('get_kalman_pll_estimates:InvalidType', 'training_scint_model must be a char or string.');
    end
    training_scint_model = char(training_scint_model); % Convert to char if needed.
    % Validate training_scint_model option.
    validModels = {'CSM','TPPSM','none'};
    if ~any(strcmpi(training_scint_model, validModels))
        error('get_kalman_pll_estimates:InvalidTrainingModel', 'training_scint_model must be one of: %s.', strjoin(validModels, ', '));
    end
    if ~isfield(kalman_pll_config, training_scint_model)
        error('get_kalman_pll_estimates:InvalidTrainingModel', 'kalman_pll_config does not have the field for training_scint_model: %s', training_scint_model);
    end

    % Validate adaptive configuration input.
    validateattributes(adaptive_config, {'struct'}, {'nonempty'}, mfilename, 'adaptive_config');
    if ~isfield(adaptive_config, 'algorithm')
        error('get_kalman_pll_estimates:MissingField', 'adaptive_config must have field "algorithm".');
    end
    if ~isfield(adaptive_config, 'hard_limited')
        error('get_kalman_pll_estimates:MissingField', 'adaptive_config must have field "hard_limited".');
    end
    % Validate the algorithm field.
    validAlgos = {'simplified','NWBP','none'};
    if ~any(strcmpi(adaptive_config.algorithm, validAlgos))
        error('get_kalman_pll_estimates:InvalidAlgorithm', 'Invalid adaptive algorithm. Allowed values: %s.', strjoin(validAlgos,', '));
    end
    if strcmpi(adaptive_config.algorithm, 'NWBP')
        error('get_kalman_pll_estimates:NotImplemented', 'The NWBP adaptive algorithm is not implemented in this version.');
    end
    % Validate the hard_limited field.
    if ~(islogical(adaptive_config.hard_limited) || isnumeric(adaptive_config.hard_limited))
        error('get_kalman_pll_estimates:InvalidHardLimited', 'adaptive_config.hard_limited must be a logical value.');
    end

    % If the simplified algorithm is chosen, require additional fields.
    if strcmpi(adaptive_config.algorithm, 'simplified')
        requiredSimplifiedFields = {'L1_C_over_N0_dBHz','sampling_interval','threshold'};
        for i = 1:length(requiredSimplifiedFields)
            if ~isfield(adaptive_config, requiredSimplifiedFields{i})
                error('get_kalman_pll_estimates:MissingField', 'adaptive_config must have field "%s" when using the simplified algorithm.', requiredSimplifiedFields{i});
            end
        end
        baseline_L1_c_over_N0_dBHz = adaptive_config.L1_C_over_N0_dBHz;
        baseline_L1_c_over_n0_linear = 10^(baseline_L1_c_over_N0_dBHz/10); % convert from dBHz to linear scale
        sampling_interval = adaptive_config.sampling_interval;
        threshold = adaptive_config.threshold;
    end

    %% Retrieve Kalman filter matrices from configuration.
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
    
    %% Preallocate output arrays.
    N = size(received_signal,1);
    state_estimates = zeros(N, numel(initial_estimates.x_hat_init));
    error_covariance_estimates = zeros(N, size(initial_estimates.P_hat_init,1), size(initial_estimates.P_hat_init,2));
    
    %% Initialize estimates.
    x_hat_project_ahead = initial_estimates.x_hat_init;
    P_hat_project_ahead = initial_estimates.P_hat_init;
    
    %% Main filtering loop.
    for step = 1:N
        if step > 1
            % Determine the measurement noise covariance adapt_R.
            % If the adaptive algorithm is 'none', no adaptation is applied.
            if strcmpi(adaptive_config.algorithm, 'none')
                adapt_R = R;
            elseif strcmpi(adaptive_config.algorithm, 'simplified')
                % Simplified adaptive update (placeholder for [3, Equation 64]).
                current_intensity = abs(received_signal(step-1,1))^2;
                estimated_L1_c_over_n0_linear = baseline_L1_c_over_n0_linear * current_intensity;
                phase_noise_variance = (1/(2 * estimated_L1_c_over_n0_linear * sampling_interval)) * ...
                                       (1 + (1/(2 * estimated_L1_c_over_n0_linear * sampling_interval)));
                adapt_R = phase_noise_variance;
                % Apply hard-limited constraint if enabled.
                if adaptive_config.hard_limited
                    if (10*log10(estimated_L1_c_over_n0_linear) < threshold)
                        adapt_R = 1e6;
                    end
                end
            else
                % Should not reach here because NWBP is not implemented.
                adapt_R = R;
            end
            
            % Compute Kalman gain using the (possibly adapted) measurement noise covariance.
            K = P_hat_project_ahead * H.' * ((H * P_hat_project_ahead * H.' + adapt_R) \ eye(size(R,1)));
            
            % Update state estimate using the phase error from the received signal.
            x_hat_update = x_hat_project_ahead + K * angle(received_signal(step-1,1) * exp(-1j * (H * x_hat_project_ahead)));
            
            % Update error covariance.
            P_hat_update = P_hat_project_ahead - K * H * P_hat_project_ahead;
        else
            x_hat_update = x_hat_project_ahead;
            P_hat_update = P_hat_project_ahead;
        end
        
        % Projection step.
        x_hat_project_ahead = F * x_hat_update + W;
        P_hat_project_ahead = F * P_hat_update * F.' + Q;
        
        % Save estimates.
        state_estimates(step,:) = x_hat_project_ahead.';
        error_covariance_estimates(step,:,:) = P_hat_project_ahead;
    end
end

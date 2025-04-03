function [state_estimates, error_covariance_estimates, L1_c_over_n0_linear_estimates] = get_kalman_pll_estimates(received_signal, kalman_pll_config, initial_estimates, training_scint_model, adaptive_config, online_mdl_learning_cfg)
% get_kalman_pll_estimates
%
% Generates Kalman filter state and error covariance estimates based on the
% received signal and the provided Kalman PLL configuration.
%
% Syntax:
%   [state_estimates, error_covariance_estimates] = get_kalman_pll_estimates(received_signal, kalman_pll_config, initial_estimates, training_scint_model, adaptive_config, online_mdl_learning_cfg)
%
% Description:
%   This function applies the Kalman filter algorithm iteratively to the
%   received signal to produce state estimates and error covariance estimates.
%   It uses the configuration settings (F, Q, H, R, W) stored in the field of
%   kalman_pll_config corresponding to the given training_scint_model.
%
%   The training_scint_model parameter can be one of:
%     'CSM'      - Use Augmentation Model trained with the Cornell 
%                  Scintillation Model (CSM).
%     'TPPSM'    - Use Augmentation Model trained with the Two Component 
%                  Power Law Phase Screen Model (TPPSM).
%     'none'     - No augmentation model).
%
%   In addition, an adaptive configuration is applied based on the provided
%   adaptive_config struct. The adaptive configuration is organized into several
%   parts:
%
%   1. Measurement Covariance Adaptation:
%      - Field: measurement_cov_adapt_algorithm
%        Allowed values: {'none', 'simplified', 'nwpr'}.
%          * 'none': No adaptive update is applied (the measurement noise covariance R is used directly).
%          * 'simplified': A simplified adaptive update is applied using the estimated carrier-to-noise ratio.
%              In this case, the field measurement_cov_adapt_algorithm_params must be provided as a struct
%              containing at least:
%                  - L1_C_over_N0_dBHz : double value representing the estimated carrier-to-noise ratio.
%          * 'nwpr': A NWPR-based update is used. When selected, measurement_cov_adapt_algorithm_params must include
%              additional fields, e.g., 'N_nwpr', 'M_nwpr' (both doubles).
%
%   2. State Covariance Adaptation:
%      - Field: states_cov_adapt_algorithm
%        Allowed values: {'none', 'matching', 'self_adaptive'}.
%          * 'none': No state covariance adaptation is applied.
%          * 'matching': An adaptive update using the covariance matching method is applied.
%              In this case, the field states_cov_adapt_algorithm_params must be provided with:
%                  - method     : string, either 'IAE' or 'RAE'.
%                  - N_nwpr: double value.
%                  - M_nwpr: double value
%          * 'self_adaptive': A self-adaptive update is applied.
%              Then, states_cov_adapt_algorithm_params must include:
%                  - init_meas_cov   : initial measurement covariance.
%                  - init_states_cov : initial states covariance.
%                  - sigma_threshold : a threshold parameter.
%
%   3. Global Sampling Interval:
%      - Field: sampling_interval (double)
%        This field is used consistently across both measurement and state adaptations.
%
%   4. Hard Limiting:
%      - Field: hard_limited (struct)
%        This struct must contain:
%                  - is_used: logical scalar indicating if hard limiting is enabled.
%                  - L1_C_over_N0_dBHz_threshold: double scalar threshold value.
%        Note: If hard limiting (hard_limited.is_used is true) is enabled, then a measurement
%        covariance adaptation algorithm must be specified (i.e., measurement_cov_adapt_algorithm must not be 'none').
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
%   adaptive_config      - A struct with adaptive configuration options for the Kalman filter.
%                          It must include the following fields:
%         .measurement_cov_adapt_algorithm
%             - Allowed values: {'none','simplified','nwpr'}.
%             - If not 'none', then the field measurement_cov_adapt_algorithm_params must be provided,
%               containing the required parameters (see above).
%         .states_cov_adapt_algorithm
%             - Allowed values: {'none','matching','self_adaptive'}.
%             - If not 'none', then the field states_cov_adapt_algorithm_params must be provided,
%               containing the required parameters (see above).
%         .sampling_interval
%             - A double specifying the global sampling interval.
%         .hard_limited
%             - A struct with the fields:
%                  is_used (logical) and L1_C_over_N0_dBHz_threshold (double).
%
%   online_mdl_learning_cfg - Struct with settings for the online model learning:
%       'is_online': Boolean flag that determines if online augmentation model learning is used.
%       For AR augmentation models ('arfit', 'aryule'):
%         'learning_method': String with one of the following methods: {'sliding_window', 'block_window', 'kalman'}.
%         'window_size' (optional): Window size for the chosen method.
%         'kalman_cfg' (optional): Kalman configuration for the online AR model parameter estimation.
%       For RBF augmentation model ('rbf'):
%         'rbf_cfg': RBF network configuration for online model learning.
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

    validate_received_signal(received_signal);
    validate_initial_estimates(initial_estimates);
    validate_training_scint_model(training_scint_model, kalman_pll_config);
    validate_adaptive_config(adaptive_config);
    validate_online_learning_cfg(online_mdl_learning_cfg, kalman_pll_config, training_scint_model);
    
    % Retrieve Kalman filter matrices from configuration.
    config_struct = kalman_pll_config.(training_scint_model);
    required_config_fields = {'F', 'Q', 'H', 'R', 'W'};
    for i = 1:length(required_config_fields)
        if ~isfield(config_struct, required_config_fields{i})
            error('get_kalman_pll_estimates:missing_config_field', 'kalman_pll_config.%s is missing the field %s.', training_scint_model, required_config_fields{i});
        end
    end
    F = config_struct.F;
    Q = config_struct.Q;
    H = config_struct.H;
    R = config_struct.R;
    W = config_struct.W;
    
    % Preallocate output arrays.
    N = size(received_signal, 1);
    state_estimates = zeros(N, numel(initial_estimates.x_hat_init));
    error_covariance_estimates = zeros(N, size(initial_estimates.P_hat_init,1), size(initial_estimates.P_hat_init,2));
    L1_c_over_n0_linear_estimates = NaN(N, 1);

    % Initialize estimates.
    x_hat_project_ahead = initial_estimates.x_hat_init;
    P_hat_project_ahead = initial_estimates.P_hat_init;
    
    if adaptive_config.hard_limited.is_used
        threshold = adaptive_config.hard_limited.L1_C_over_N0_dBHz_threshold;
    end

    % If the simplified adaptive algorithm is chosen, compute baseline parameters.
    if strcmpi(adaptive_config.measurement_cov_adapt_algorithm, 'simplified')
        baseline_L1_C_over_N0_dBHz = adaptive_config.measurement_cov_adapt_algorithm_params.L1_C_over_N0_dBHz;
        baseline_L1_c_over_n0_linear = 10^(baseline_L1_C_over_N0_dBHz / 10);
        sampling_interval = adaptive_config.sampling_interval;
    end
    if strcmpi(adaptive_config.measurement_cov_adapt_algorithm, 'nwpr')
        N_nwpr = adaptive_config.measurement_cov_adapt_algorithm_params.N_nwpr;
        M_nwpr = adaptive_config.measurement_cov_adapt_algorithm_params.M_nwpr;
        sampling_interval = adaptive_config.sampling_interval;
        z_vec = received_signal(1:M_nwpr);

        % Compute NP for the initial window:
        NBP = (sum(real(z_vec)))^2 + (sum(imag(z_vec)))^2;
        WBP = sum(abs(z_vec).^2);
        current_NP = NBP / WBP;
        vec_size = N_nwpr / M_nwpr;
        % Initialize NP_vec by filling it with the computed current_NP value:
        NP_vec = repmat(current_NP, vec_size, 1);
        NP_array = zeros(N,1);
    end
    
    %% Main filtering loop.
    for step = 1:N
        if step > 1
            % Determine the measurement noise covariance adapt_R.
            if strcmpi(adaptive_config.measurement_cov_adapt_algorithm, 'none')
                adapt_R = R;
            elseif strcmpi(adaptive_config.measurement_cov_adapt_algorithm, 'simplified')
                current_intensity = abs(received_signal(step-1,1))^2;
                estimated_L1_c_over_n0_linear = current_intensity * baseline_L1_c_over_n0_linear;
                phase_noise_variance = (1/(2 * estimated_L1_c_over_n0_linear * sampling_interval)) * ...
                                       (1 + (1/(2 * estimated_L1_c_over_n0_linear * sampling_interval)));
                adapt_R = phase_noise_variance;
                
                % Store the estimate
                L1_c_over_n0_linear_estimates(step,1) = estimated_L1_c_over_n0_linear;
            elseif strcmpi(adaptive_config.measurement_cov_adapt_algorithm, 'nwpr')
                NBP = (sum(real(z_vec)))^2 + (sum(imag(z_vec)))^2;
                WBP = sum(abs(z_vec).^2);

                NP = NBP / WBP;
                NP_vec = [NP; NP_vec(1:end-1)];
                mu_hat = (M_nwpr / N_nwpr)  * sum(NP_vec);
                estimated_L1_c_over_n0_linear = (1/sampling_interval) * (mu_hat) / (M_nwpr - mu_hat);

                phase_noise_variance = (1/(2 * estimated_L1_c_over_n0_linear * sampling_interval)) * ...
                                       (1 + (1/(2 * estimated_L1_c_over_n0_linear * sampling_interval)));
                adapt_R = phase_noise_variance;

                % Update
                z_vec = [received_signal(step-1,1) * exp(-1j * H * x_hat_project_ahead); z_vec(1:end-1)];
                NP_array(step,1) = NP;
                L1_c_over_n0_linear_estimates(step,1) = estimated_L1_c_over_n0_linear;
            else
                adapt_R = R;
            end
            % Apply hard-limited constraint if enabled.
            if adaptive_config.hard_limited.is_used
                if (10 * log10(estimated_L1_c_over_n0_linear) < threshold)
                    adapt_R = 1e6;
                end
            end
            if online_mdl_learning_cfg.is_online
                [F, Q, W] = update_filter_matrices(F, Q, W, step, state_estimates, online_mdl_learning_cfg, kalman_pll_config.(training_scint_model).augmentation_model_initializer);
            end
            % Compute Kalman gain.
            K = P_hat_project_ahead * H.' * ((H * P_hat_project_ahead * H.' + adapt_R) \ eye(size(R, 1)));
            
            % Update state estimate.
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
        state_estimates(step, :) = x_hat_project_ahead.';
        error_covariance_estimates(step, :, :) = P_hat_project_ahead;
    end
end

%% Helper Functions

function validate_received_signal(received_signal)
    validateattributes(received_signal, {'numeric'}, {'nonempty', '2d'}, mfilename, 'received_signal');
end

function validate_initial_estimates(initial_estimates)
    validateattributes(initial_estimates, {'struct'}, {'nonempty'}, mfilename, 'initial_estimates');
    if ~isfield(initial_estimates, 'x_hat_init') || ~isfield(initial_estimates, 'P_hat_init')
        error('get_kalman_pll_estimates:missing_field', 'initial_estimates must have fields "x_hat_init" and "P_hat_init".');
    end
    validateattributes(initial_estimates.x_hat_init, {'numeric'}, {'nonempty', 'vector'}, mfilename, 'initial_estimates.x_hat_init');
    n = numel(initial_estimates.x_hat_init);
    validateattributes(initial_estimates.P_hat_init, {'numeric'}, {'nonempty', '2d', 'size', [n n]}, mfilename, 'initial_estimates.P_hat_init');
end

function validate_training_scint_model(training_scint_model, kalman_pll_config)
    if ~(ischar(training_scint_model) || isstring(training_scint_model))
        error('get_kalman_pll_estimates:invalid_type', 'training_scint_model must be a char or string.');
    end
    training_scint_model = char(training_scint_model);
    valid_models = {'CSM', 'TPPSM', 'none'};
    if ~any(strcmpi(training_scint_model, valid_models))
        error('get_kalman_pll_estimates:invalid_training_model', 'training_scint_model must be one of: %s.', strjoin(valid_models, ', '));
    end
    if ~isfield(kalman_pll_config, training_scint_model)
        error('get_kalman_pll_estimates:invalid_training_model', 'kalman_pll_config does not have the field for training_scint_model: %s', training_scint_model);
    end
end

function validate_adaptive_config(adaptive_config)
    validateattributes(adaptive_config, {'struct'}, {'nonempty'}, mfilename, 'adaptive_config');

    % Required main fields
    required_main_fields = {'measurement_cov_adapt_algorithm', 'states_cov_adapt_algorithm', 'sampling_interval', 'hard_limited'};
    for i = 1:numel(required_main_fields)
        if ~isfield(adaptive_config, required_main_fields{i})
            error('validate_adaptive_config:missing_field', ...
                  'adaptive_config must contain the field "%s".', required_main_fields{i});
        end
    end

    % Validate sampling_interval
    validateattributes(adaptive_config.sampling_interval, {'numeric'}, {'scalar', 'positive'}, ...
                       mfilename, 'adaptive_config.sampling_interval');

    % Validate hard_limited
    hl = adaptive_config.hard_limited;
    validateattributes(hl, {'struct'}, {'nonempty'}, mfilename, 'adaptive_config.hard_limited');
    if ~isfield(hl,'is_used')
        error('validate_adaptive_config:missing_hard_limited_field', ...
                  'hard_limited must have field "%s".', 'is_used');
    end
    validateattributes(hl.is_used, {'logical', 'numeric'}, {'scalar'}, mfilename, 'hard_limited.is_used');
    if hl.is_used
        if ~isfield(hl, 'L1_C_over_N0_dBHz_threshold')
            error('validate_adaptive_config:missing_hard_limited_field', ...
                  'hard_limited must have field "%s".', 'L1_C_over_N0_dBHz_threshold');
        end
        validateattributes(hl.L1_C_over_N0_dBHz_threshold, {'numeric'}, {'scalar'}, ...
                       mfilename, 'hard_limited.L1_C_over_N0_dBHz_threshold');
    end
    
    % Validate measurement_cov_adapt_algorithm
    valid_meas_algos = {'none', 'simplified', 'nwpr'};
    algo = adaptive_config.measurement_cov_adapt_algorithm;
    assert(any(strcmpi(algo, valid_meas_algos)), ...
        'validate_adaptive_config:invalid_algorithm', ...
        'Invalid measurement_cov_adapt_algorithm. Allowed values are: %s', strjoin(valid_meas_algos, ', '));

    if hl.is_used && strcmpi(algo, 'none')
        error('validate_adaptive_config:hard_limited_without_adaptive', ...
              ['The Hard limited option should only be configured to be used ' ...
              'when a measurement covariance algorithm is also used.'])
    end
    if ~strcmpi(algo, 'none')
        if ~isfield(adaptive_config, 'measurement_cov_adapt_algorithm_params')
            error('validate_adaptive_config:missing_field', ...
                  'Field "measurement_cov_adapt_algorithm_params" is required for algorithm "%s".', algo);
        end
        p = adaptive_config.measurement_cov_adapt_algorithm_params;
        validateattributes(p, {'struct'}, {'nonempty'}, mfilename, 'measurement_cov_adapt_algorithm_params');
        % Shared required fields
        base_fields = {'L1_C_over_N0_dBHz'};
        if strcmpi(algo, 'simplified')
            check_fields(p, base_fields, 'measurement_cov_adapt_algorithm_params');
        elseif strcmpi(algo, 'nwpr')
            extra_fields = {'N_nwpr','M_nwpr'};
            check_fields(p, [base_fields, extra_fields], 'measurement_cov_adapt_algorithm_params');
        end
    end

    % Validate states_cov_adapt_algorithm
    valid_states_algos = {'none', 'matching', 'self_adaptive'};
    s_algo = adaptive_config.states_cov_adapt_algorithm;
    assert(any(strcmpi(s_algo, valid_states_algos)), ...
        'validate_adaptive_config:invalid_states_algorithm', ...
        'Invalid states_cov_adapt_algorithm. Allowed values are: %s', strjoin(valid_states_algos, ', '));

    if strcmpi(s_algo, 'matching')
        if ~isfield(adaptive_config, 'states_cov_adapt_algorithm_params')
            error('validate_adaptive_config:missing_field', ...
                  'Field "states_cov_adapt_algorithm_params" is required for algorithm "%s".', s_algo);
        end
        p = adaptive_config.states_cov_adapt_algorithm_params;
        required_fields = {'method', 'window_size'};
        check_fields(p, required_fields, 'states_cov_adapt_algorithm_params');
        if ~any(strcmpi(p.method, {'IAE', 'RAE'}))
            error('validate_adaptive_config:invalid_method', ...
                  'Invalid method in states_cov_adapt_algorithm_params. Must be either "IAE" or "RAE".');
        end
    elseif strcmpi(s_algo, 'self_adaptive')
        required_fields = {'init_meas_cov', 'init_states_cov', 'sigma_threshold'};
        p = adaptive_config.states_cov_adapt_algorithm_params;
        check_fields(p, required_fields, 'states_cov_adapt_algorithm_params');
    end
end

function check_fields(structure, field_list, context)
    for i = 1:numel(field_list)
        if ~isfield(structure, field_list{i})
            error('validate_adaptive_config:missing_field', ...
                  '%s must contain the field "%s".', context, field_list{i});
        end
    end
end

function validate_online_learning_cfg(online_mdl_learning_cfg, kalman_pll_config, training_scint_model)
    validateattributes(online_mdl_learning_cfg, {'struct'}, {'nonempty'}, mfilename, 'online_mdl_learning_cfg');
    
    if ~isfield(online_mdl_learning_cfg, 'is_online')
        error('get_kalman_pll_estimates:missing_field', ...
              'online_mdl_learning_cfg must have a field "is_online" that indicates if online model learning is enabled.');
    end
    validateattributes(online_mdl_learning_cfg.is_online, {'logical'}, {}, mfilename, 'online_mdl_learning_cfg.is_online');
    
    if online_mdl_learning_cfg.is_online
        % Get the offline augmentation model identifier.
        if ~isfield(kalman_pll_config.(training_scint_model), 'augmentation_model_initializer')
            error('get_kalman_pll_estimates:missing_field', ...
                  'kalman_pll_config must have a field "augmentation_model_initializer" specifying the offline augmentation model.');
        end
        if ~isfield(kalman_pll_config.(training_scint_model).augmentation_model_initializer, 'id')
            error('get_kalman_pll_estimates:missing_field', ...
                  'kalman_pll_config.augmentation_model_initializer must contain the field "id" that identifies the augmentation model used in offline training.');
        end
        offline_id = kalman_pll_config.(training_scint_model).augmentation_model_initializer.id;
        if ~(ischar(offline_id) || isstring(offline_id))
            error('get_kalman_pll_estimates:invalid_type', ...
                  'The field kalman_pll_config.augmentation_model_initializer.id must be a char or string.');
        end
        offline_id = char(offline_id);
        
        switch lower(offline_id)
            case {'arfit', 'aryule'}
                % For AR models, require the fields directly in online_mdl_learning_cfg.
                if ~isfield(online_mdl_learning_cfg, 'learning_method')
                    error('get_kalman_pll_estimates:missing_field', ...
                          'For offline augmentation model id "%s", online_mdl_learning_cfg must have a field "learning_method".', offline_id);
                end
                validateattributes(online_mdl_learning_cfg.learning_method, {'char','string'}, {'nonempty'}, mfilename, 'online_mdl_learning_cfg.learning_method');
                method = lower(char(online_mdl_learning_cfg.learning_method));
                if ismember(method, {'sliding_window', 'block_window'})
                    if ~isfield(online_mdl_learning_cfg, 'window_size')
                        error('get_kalman_pll_estimates:missing_field', ...
                              'For learning method "%s" with offline model id "%s", online_mdl_learning_cfg must include a "window_size" field.', method, offline_id);
                    end
                    validateattributes(online_mdl_learning_cfg.window_size, {'numeric'}, {'scalar','positive'}, mfilename, 'online_mdl_learning_cfg.window_size');
                end
                if ismember(method, 'kalman')
                    error('get_kalman_pll_estimates:NotSupportedLearningMethod', 'Adapting the AR coefficients using a Kalman filter is not yet supported.')
                end
            case 'rbf'
                error('get_kalman_pll_estimates:NotSupported', 'RBF online learning is not yet supported.');
            case 'kinematic'
                error('Online model estimation for ''second_wiener_mdl'' is not yet supported.');
            case 'none'
                % Does nothing
            otherwise
                error('get_kalman_pll_estimates:invalid_augmentation_id', ...
                      'Unsupported offline augmentation model id "%s" found in kalman_pll_config.augmentation_model_initializer.id', offline_id);
        end
    end
end


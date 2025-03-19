function [phase_estimates, frequency_estimates] = get_fll_estimates(fll_config, received_signal)
% GET_FLL_ESTIMATES
%
% Syntax:
%   [phase_estimates, frequency_estimates] = get_fll_estimates(fll_config, received_signal)
%
% Description:
%   Implements a Frequency Locked Loop (FLL) to estimate the phase and
%   frequency of a received complex signal. The function processes the signal
%   sample-by-sample based on the previous phase estimate, a phase error 
%   discriminator, and a digital loop filter to update the estimates.
%
%   The digital loop filter design used in this function is derived from 
%   the approach presented in [1]. In particular, the manual scaling parameter 
%   kd (provided in fll_config.loop_filter_cfg.kd) is fully characterized in 
%   [1] and depends on both the selected discriminator and the signal-to-noise 
%   ratio (SNR). The SNR can be computed from the carrier-to-noise ratio by 
%   dividing its linear form by the receiver bandwidth, which is assumed to be 
%   half of the receiver's sampling frequency.
%
% Inputs:
%   fll_config       - Structure with FLL configuration parameters including:
%                      discriminator_type, sampling_interval,
%                      initial_frequency_estimate_boundaries, and loop_filter_cfg.
%   received_signal  - Column vector of complex samples representing the signal.
%
% Outputs:
%   phase_estimates     - Column vector of estimated phase values (radians).
%   frequency_estimates - Column vector of estimated frequency values (rad/s).
%
% Notes:
%   - Currently only 'atan2' discriminator is implemented; other types trigger
%     errors. You may add more later.
%   - We assume here that 'kd' is externally provided. Future versions could
%     compute kd based on SNR.
%
% Example:
%   % Example usage:
%   fll_config.discriminator_type = 'atan2';
%   fll_config.sampling_interval = 1e-3;
%   fll_config.initial_frequency_estimate_boundaries = [-50, 50];
%   fll_config.loop_filter_cfg.order = 2;
%   fll_config.loop_filter_cfg.noise_bandwidth = 5;
%   fll_config.loop_filter_cfg.damping_mode = 'critically_damped';
%   fll_config.loop_filter_cfg.kd = 0.8;
%   t = (0:1e-3:1)';
%   received_signal = exp(1j*(2*pi*10*t + 0.1*randn(size(t))));
%   [phase_est, freq_est] = get_fll_estimates(fll_config, received_signal);
%
% Dependencies:
%   design_filter, validate_fll_config.
%
% References:
% [1] Curran, James T. & Lachapelle, Gérard & Murphy, Colin. (2012). 
%     Improving the Design of Frequency Lock Loops for GNSS Receivers.
%     IEEE Transactions on Aerospace and Electronic Systems. 
%     doi: 10.1109/TAES.2012.6129674
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

   
    % Validate FLL configuration
    validate_fll_config(fll_config);

    % Extract FLL parameters
    sampling_interval = fll_config.sampling_interval;
    discriminator_type = char(fll_config.discriminator_type);
    init_freq_est_boundaries = fll_config.initial_frequency_estimate_boundaries;
    
    % Choose an initial frequency estimate in Hz
    init_freq_est_Hz = unifrnd(init_freq_est_boundaries(1), init_freq_est_boundaries(2));
    initial_phase_rad = 0;

    % Design loop filter and obtain coefficients
    [filter_den, filter_num] = design_filter(fll_config);

    % Initialize filter state vector
    filter_state_length = fll_config.loop_filter_cfg.order - 1;
    filter_state = zeros(filter_state_length, 1);

    % Number of samples
    N = length(received_signal);

    % Pre-allocate outputs
    phase_estimates = zeros(N, 1);
    frequency_estimates = zeros(N, 1);

    % Convert initial frequency from Hz to rad/s
    init_freq_est_rad_s = 2 * pi * init_freq_est_Hz;

    % Set initial conditions at sample 1
    phase_estimates(1) = initial_phase_rad + 2*pi*init_freq_est_rad_s * sampling_interval;
    frequency_estimates(1) = init_freq_est_rad_s;

    % Initialize the previous error signal using sample 1
    error_signal_prev = received_signal(1) * exp(-1j * phase_estimates(1));

    % Loop from sample 2 to the end
    for step = 2:N
        % De-rotate the current received sample using the previous phase estimate
        error_signal = received_signal(step) * exp(-1j * phase_estimates(step-1));

        % Compute the discriminator output (phase_error in radians)
        switch discriminator_type
            case 'atan2'
                dot_val = real(error_signal) * real(error_signal_prev) + ...
                          imag(error_signal) * imag(error_signal_prev);
                cross_val = imag(error_signal) * real(error_signal_prev) - ...
                            real(error_signal) * imag(error_signal_prev);
                phase_error = atan2(dot_val, cross_val);
            otherwise
                error('Unsupported or unimplemented discriminator type: %s', discriminator_type);
        end

        % Convert phase error to a frequency error (rad/s)
        freq_error = phase_error / sampling_interval;

        % Filter the frequency error to obtain the frequency correction
        [freq_correction, filter_state] = filter(filter_num, filter_den, freq_error, filter_state);

        % Update the frequency estimate (integrate the frequency correction)
        frequency_estimates(step) = frequency_estimates(step-1) + freq_correction;

        % Integrate the frequency estimate to update phase
        phase_estimates(step) = phase_estimates(step-1) + sampling_interval * frequency_estimates(step);

        % Update the previous error signal for the next iteration
        error_signal_prev = error_signal;
    end
end

function [filter_den, filter_num] = design_filter(fll_config)
% design_filter returns digital loop filter coefficients based on the
% loop filter configuration in loop_cfg.
%
% The filter coefficients are determined by the filter order and, for a 
% second order filter, by the damping mode. The manual scaling parameter kd
% (located in loop_cfg.kd) is used in the design.
%
% Reference
% [1] Curran, James T. & Lachapelle, Gérard & Murphy, Colin. (2012). 
%     Improving the Design of Frequency Lock Loops for GNSS Receivers.
%     IEEE Transactions on Aerospace and Electronic Systems. 
%     doi: 10.1109/TAES.2012.6129674
    
    order = fll_config.loop_filter_cfg.order;
    kd = fll_config.loop_filter_cfg.kd;  % manual input parameter for scaling
    bw = fll_config.loop_filter_cfg.noise_bandwidth;
    tl = fll_config.sampling_interval;
    switch order
        case 1
            % First order loop filter: no damping mode selection is needed.
            % Use the manual parameter kd directly.
            A0 = (2*bw)/(kd*(1+bw));
            filter_num = A0;
            filter_den = 1;
        case 2
            % Second order loop filter: design depends on the damping mode.
            if strcmpi(fll_config.loop_filter_cfg.damping_mode, 'critically_damped')
                % Equation 25 of [1]
                A0 = (1 - exp((1/6) * bw * tl * (2 * bw * tl - 3*pi))) / (kd * tl); 
                % Equation 26 of [1]
                A1 = (exp(-(1/2) * bw * tl * pi) * (exp((1/6)* (bw * tl)^2) - exp((1/4) * bw * tl * pi))^2) / (kd * tl^2);
            elseif strcmpi(fll_config.loop_filter_cfg.damping_mode, 'underdamped')
                % Equation 29 of [1]
                gamma = (1/5) * bw * tl * (bw * tl - (41/(4*pi)));
                % Equation 27 of [1]
                A0 = (1 - exp(2*gamma)) / (kd * tl);
                % Equation 28 of [1]
                A1 = (1 + exp(2*gamma) - 2 * exp(gamma) * cos(gamma)) / (kd * tl^2);
            else
                error('Unsupported damping mode: %s', fll_config.loop_filter_cfg.damping_mode);
            end
            filter_num = [A0 + A1 * tl, -A0];
            filter_den = [1, -1];
        otherwise
            error('Unsupported loop filter order: %d', order);
    end
end

function validate_fll_config(fll_config)
    % Validate the overall configuration struct.
    validateattributes(fll_config, {'struct'}, {'nonempty'}, mfilename, 'fll_config');
    
    %% Validate fll_config fields
    
    % 1. discriminator_type: must be a nonempty char or string and one of the allowed values.
    requiredField = 'discriminator_type';
    if ~isfield(fll_config, requiredField)
        error('fll_config must include the field "%s".', requiredField);
    end
    validateattributes(fll_config.discriminator_type, {'char', 'string'}, {'nonempty'}, mfilename, requiredField);
    valid_discriminators = {'atan2', 'cp', 'atan', 'ddcp'};
    if ~ismember(char(fll_config.discriminator_type), valid_discriminators)
        error('fll_config.discriminator_type must be one of: %s', strjoin(valid_discriminators, ', '));
    end
    
    % 2. sampling_interval: numeric scalar and positive.
    requiredField = 'sampling_interval';
    if ~isfield(fll_config, requiredField)
        error('fll_config must include the field "%s".', requiredField);
    end
    validateattributes(fll_config.sampling_interval, {'numeric'}, {'scalar', 'positive', 'real'}, mfilename, requiredField);
    
    % 3. initial_frequency_estimate_boundaries: numeric vector with two elements.
    requiredField = 'initial_frequency_estimate_boundaries';
    if ~isfield(fll_config, requiredField)
        error('fll_config must include the field "%s".', requiredField);
    end
    validateattributes(fll_config.initial_frequency_estimate_boundaries, {'numeric'}, ...
        {'vector', 'numel', 2, 'real'}, mfilename, requiredField);
    bounds = fll_config.initial_frequency_estimate_boundaries;
    if bounds(1) >= bounds(2)
        error('The lower bound of initial_frequency_estimate_boundaries must be less than the upper bound.');
    end
    
    % 4. loop_filter_cfg: must be a nonempty struct.
    requiredField = 'loop_filter_cfg';
    if ~isfield(fll_config, requiredField)
        error('fll_config must include the field "%s".', requiredField);
    end
    validateattributes(fll_config.loop_filter_cfg, {'struct'}, {'nonempty'}, mfilename, requiredField);
    loop_cfg = fll_config.loop_filter_cfg;
    
    %%% Validate loop_filter_cfg fields
    
    % a) order: numeric scalar, integer and either 1 or 2.
    subField = 'order';
    if ~isfield(loop_cfg, subField)
        error('loop_filter_cfg must include the field "%s".', subField);
    end
    validateattributes(loop_cfg.order, {'numeric'}, {'scalar', 'integer'}, mfilename, ['loop_filter_cfg.' subField]);
    if ~ismember(loop_cfg.order, [1, 2])
        error('loop_filter_cfg.order must be either 1 or 2.');
    end
    
    % b) noise_bandwidth: numeric scalar, positive.
    subField = 'noise_bandwidth';
    if ~isfield(loop_cfg, subField)
        error('loop_filter_cfg must include the field "%s".', subField);
    end
    validateattributes(loop_cfg.noise_bandwidth, {'numeric'}, {'scalar', 'positive', 'real'}, mfilename, ['loop_filter_cfg.' subField]);
    
    % c) kd: manual scaling parameter; numeric scalar, positive.
    subField = 'kd';
    if ~isfield(loop_cfg, subField)
        error('loop_filter_cfg must include the field "%s".', subField);
    end
    validateattributes(loop_cfg.kd, {'numeric'}, {'scalar', 'positive', 'real'}, mfilename, ['loop_filter_cfg.' subField]);
    
    % d) If order is 2, validate damping_mode: nonempty char or string from allowed values.
    if loop_cfg.order == 2
        subField = 'damping_mode';
        if ~isfield(loop_cfg, subField)
            error('For a second order filter, loop_filter_cfg must include the field "%s".', subField);
        end
        validateattributes(loop_cfg.damping_mode, {'char', 'string'}, {'nonempty'}, mfilename, ['loop_filter_cfg.' subField]);
        valid_damping_modes = {'critically_damped', 'underdamped'};
        if ~ismember(char(loop_cfg.damping_mode), valid_damping_modes)
            error('loop_filter_cfg.damping_mode must be one of: %s', strjoin(valid_damping_modes, ', '));
        end
    end
end

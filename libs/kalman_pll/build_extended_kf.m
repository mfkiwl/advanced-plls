function [F, Q, Hj_handle, R, W] = build_extended_kf(F_los, Q_los, aug_data, general_config, sampling_interval)
% build_extended_kf
%
% Constructs the full Kalman Filter matrices and returns a handle to the 
% measurement Jacobian function.
%
% Inputs:
%   F_los            - LOS dynamics state transition matrix.
%   Q_los            - LOS dynamics process noise covariance.
%   aug_data         - Struct returned by compute_augmentation containing the 
%                      augmentation type.
%   general_config   - Configuration struct containing, among others,
%                      C_over_N0_array_dBHz.
%   sampling_interval- The sampling interval.
%
% Outputs:
%   F, Q, R, W       - The state transition, process noise, measurement noise,
%                      and additional KF matrices.
%   Hj_handle        - A function handle to compute the measurement Jacobian.
%
% Notes:
%   - In this version, the amplitude is not part of the KF state but is 
%     provided by an external estimator.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Default measurement noise covariance: Adjust using get_phase_variances.
    variances = get_phase_variances(general_config.C_over_N0_array_dBHz, sampling_interval);
    I_Q_separate_variances = repelem(variances, 2) / 2;
    R = I_Q_separate_variances;

    % Choose matrices and measurement Jacobian based on augmentation type.
    switch lower(string(aug_data.augmentation_type))
        case 'none'
            % No augmentation is applied.
            F = F_los;
            Q = Q_los;
            states_amount = size(F_los, 1);  % Number of KF states (excluding amplitude)
            W = zeros(states_amount, 1);

            % Create a parameter structure for the measurement Jacobian function.
            % For the 'none' augmentation, assume that the KF state vector x contains
            % (for example) phase as its first element. (If your state includes more,
            % assign the proper index and extend the Jacobian accordingly.)
            meas_param.idx.phase = 1;      % x(1) is assumed to be phase.
            meas_param.states_amount = states_amount;
            
            % Create a handle that binds the meas_param structure. The measurement
            % Jacobian will require two inputs: the external amplitude and the KF state.
            Hj_handle = @(external_amplitude, x) measurement_jacobian(external_amplitude, x, meas_param);

        case {'arfit', 'aryule', 'kinematic', 'arima'}
            % Future augmentation types to be implemented.
            error('MATLAB:NotImplementedAugmentationModel', ...
                'KF augmentation type %s is not supported in extended KF builder yet.', aug_data.augmentation_type);

        otherwise
            error('MATLAB:InvalidAugmentationModel', ...
                'KF augmentation type %s is not supported in extended KF builder.', aug_data.augmentation_type);
    end

end

%% Helper Function: Measurement Jacobian
function Hj = measurement_jacobian(external_amplitude, x, param)
% measurement_jacobian Computes the Jacobian matrix for the measurement function.
%
% In this implementation, the measurement model is assumed to be
%
%      h(x) = [ I; Q ]  = [ external_amplitude*cos(phase);
%                           external_amplitude*sin(phase) ]
%
% where *phase* is a state variable (indexed by param.idx.phase) in the KF state x,
% and the amplitude is provided externally.
%
% Inputs:
%   external_amplitude - The amplitude estimate provided by an external estimator.
%   x                  - The KF state vector.
%   param              - A struct with fields:
%                          .idx.phase     - The index of the phase within x.
%                          .states_amount - The total number of states.
%
% Output:
%   Hj    - The measurement Jacobian matrix.
%
% For this model, note that the measurement does not depend on amplitude (since
% it's externally provided), so the derivatives with respect to the KF state are:
%     dI/dphase = -external_amplitude * sin(phase)
%     dQ/dphase =  external_amplitude * cos(phase)
%
% The Jacobian is two rows (I and Q) by the number of KF states.
    
    % Extract the phase from the state vector.
    phase_estimate = x(param.idx.phase);
    states_amount  = param.states_amount;
    
    % Compute partial derivatives with respect to the phase.
    dI_dphase = -external_amplitude * sin(phase_estimate);
    dQ_dphase =  external_amplitude * cos(phase_estimate);
    
    % Initialize the Jacobian as a 2 x states_amount matrix.
    % Here, only the column corresponding to phase (param.idx.phase) is non-zero.
    Hj = zeros(2, states_amount);
    Hj(:, param.idx.phase) = [dI_dphase; dQ_dphase];
    
    % Future extensions: If more state variables affect the measurement, compute 
    % and assign their derivatives here.
end

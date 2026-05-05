function [F, Q, R, W] = build_unscented_kf(F_los, Q_los, aug_data, general_config, sampling_interval)
% build_unscented_kf
%
% Constructs the full unscented Kalman filter matrices.
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
%                      and additive state-transition vector.
% Notes:
%   - In this version, the amplitude is not part of the KF state but is 
%     provided by an external estimator.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Default measurement noise covariance: Adjust using get_phase_variances.
    variances = 1/(10.^(general_config.C_over_N0_array_dBHz ./ 10) .* sampling_interval);
    I_Q_separate_variances = repelem(variances, 2) / 2;
    R = diag(I_Q_separate_variances);

    % Choose matrices and measurement Jacobian based on augmentation type.
    switch lower(string(aug_data.augmentation_type))
        case 'none'
            % No augmentation is applied.
            F = F_los;
            Q = Q_los;
            W = zeros(size(F_los, 1), 1);

        case 'arfit-complex-field'
            F_aug = aug_data.F_aug;
            Q_aug = aug_data.Q_aug;
            intercept_vector = aug_data.intercept;
            states_amount = aug_data.augInfo.states_amount;
            model_order = aug_data.augInfo.model_order;

            if states_amount ~= 2
                error('MATLAB:InvalidAugmentationModel', ...
                    'The arfit-complex-field augmentation model must contain exactly real and imaginary VAR states.');
            end

            F = blkdiag(F_los, F_aug);
            Q = blkdiag(Q_los, Q_aug);
            W = [zeros(size(F_los, 1), 1); intercept_vector; zeros(states_amount * (model_order - 1), 1)];

        case {'arfit', 'aryule', 'kinematic', 'arima'}
            % Future augmentation types to be implemented.
            error('MATLAB:NotImplementedAugmentationModel', ...
                'KF augmentation type %s is not supported in unscented KF builder yet.', aug_data.augmentation_type);

        otherwise
            error('MATLAB:InvalidAugmentationModel', ...
                'KF augmentation type %s is not supported in unscented KF builder.', aug_data.augmentation_type);
    end

end

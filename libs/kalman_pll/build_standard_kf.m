function [F, Q, H, R, W] = build_standard_kf(F_los, Q_los, aug_data, general_config, sampling_interval)
% build_standard_kf constructs the full KF matrices for the standard variant.
%
% Inputs:
%   F_los, Q_los  - LOS dynamics matrices.
%   aug_data       - Struct returned by compute_augmentation.
%   general_config- Configuration struct containing, among others, C_over_N0_array_dBHz.
%   sampling_interval - The sampling interval.
%
% Outputs:
%   F, Q, H, R, W - The state transition, process noise, measurement, measurement noise,
%                   and additional matrices for the Kalman filter.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

    % Default measurement noise covariance:
    R = diag(compute_phase_variances(general_config.C_over_N0_array_dBHz, sampling_interval));

    switch lower(string(aug_data.augmentation_type))
        case 'none'
            % No augmentation is applied.
            F = F_los;
            Q = Q_los;
            % Construct H for the LOS part.
            H = [1, zeros(1, size(F_los,1)-1)];
            W = zeros(size(F_los,1), 1);

        case {'arfit', 'aryule'}
            % For ARFIT or ARYULE, we assume the existence of a helper function
            % construct_kalman_matrices that fuses the LOS matrices with the VAR–based augmentation.
            % The function is assumed to have the following signature:
            %
            % [F, Q, H, R, W] = construct_kalman_matrices(F_los, Q_los, F_aug, Q_aug, intercept, states_amount, model_order, C_over_N0_array_dBHz, sampling_interval)
            %
            F_aug = aug_data.F_aug;
            Q_aug = aug_data.Q_aug;
            % Extract extra information (for ARYULE, intercept is set to 0).
            intercept = aug_data.intercept;
            states_amount = aug_data.augInfo.states_amount;
            model_order = aug_data.augInfo.model_order;
            [F, Q, H, R, W] = construct_kalman_matrices( ...
                F_los, Q_los, F_aug, Q_aug, intercept, states_amount, model_order, ...
                general_config.C_over_N0_array_dBHz, sampling_interval);

        case 'kinematic'
            % For the kinematic augmentation, simply concatenate the LOS and Wiener augmentation matrices.
            F = blkdiag(F_los, aug_data.F_aug);
            Q = blkdiag(Q_los, aug_data.Q_aug);
            % Construct H by selecting the first element of the LOS part and the first element of the augmentation part.
            n_los = size(F_los, 1);
            n_aug = size(aug_data.F_aug, 1);
            H = [1, zeros(1, n_los-1), 1, zeros(1, n_aug-1)];
            W = zeros(n_los + n_aug, 1);

        case 'arima'
            % For the ARIMA model, combine the LOS and ARIMA matrices.
            F = blkdiag(F_los, aug_data.F_aug);
            Q = blkdiag(Q_los, aug_data.Q_aug);
            % Use the H vector computed in the ARIMA augmentation branch.
            H = [1, zeros(1, size(F_los,1)-1), aug_data.H_aug];
            W = zeros(size(F,1), 1);

        otherwise
            error('MATLAB:UndefinedKFType', 'KF augmentation type %s is not supported in standard KF builder.', aug_data.augmentation_type);
    end
end

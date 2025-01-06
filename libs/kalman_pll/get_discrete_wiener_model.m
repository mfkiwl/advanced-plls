function [state_transition_matrix, covariance_matrix] = get_discrete_wiener_model(L, M, sampling_interval, ...
    sigma_array, delta_array)
% get_discrete_wiener_model
% Generates the state transition matrix and noise covariance matrix for
% line-of-sight (LOS) dynamics in a multi-frequency system.
%
% Syntax:
%   [state_transition_matrix, covariance_matrix] = get_discrete_wiener_model(L, M, sampling_interval, sigma, delta)
%
% Description:
%   This function computes the state transition matrix (state_transition_matrix) and the noise
%   covariance matrix (covariance_matrix) for a generic multi-frequency LOS dynamics model
%   based on the mathematical framework described in [Section 2, 1],
%   which has been inspired in the mathematical development prosed in 
%   [Section 2.4, 2].
%
% Inputs:
%   L                 - Number of carriers included in the state-space model.
%   M                 - Order of the Wiener process, indicating the degree
%                       of dynamics considered in the model.
%   sampling_interval - Time interval between samples (in seconds).
%   sigma_array       - Vector containing variances for the noise covariance
%                       matrix. Each variance reflects uncertainty related
%                       to clock drift and LOS dynamics for both satellite
%                       and receiver.
%   delta_array         - Vector of normalized frequencies for each carrier
%                       relative to the reference carrier. For example:
%                       delta = [f_L1/f_L1, f_L2/f_L1, f_L5/f_L1] for 
%                       multi-frequency tracking, or delta = 1 for 
%                       single-frequency tracking.
%
% Outputs:
%   state_transition_matrix                - State transition matrix of size (M+L-1) x (M+L-1),
%                       representing the dynamics of the state over time.
%   covariance_matrix                - Noise covariance matrix of size (M+L-1) x (M+L-1),
%                       capturing the uncertainty in the system's states.
%
% Notes:
%   - state_transition_matrix is derived using the matrix exponential of 
%     the system dynamics matrix (Fw) scaled by the sampling interval.
%   - covariance_matrix is computed using the integral of the state 
%     transition matrix over time, incorporating the Wiener process noise
%     covariance.
%   - The reference carrier frequency used in `delta` normalization must 
%     match the system's chosen base frequency (e.g., L1).
%
% Example:
% % Generate state_transition_matrix and covariance_matrix for a 3-carrier system with:
% % - Second-order Wiener process
% % - Sampling interval of 0.01 seconds
% % - Variances: [0, 0, 0, 0, 1e-2]
% % - Frequencies normalized as: delta = [1, 0.9, 0.8]
% L = 3;             % Number of carriers
% M = 3;             % Order of the Wiener process
% sampling_interval = 0.01; % Sampling interval in seconds
% sigma_array = [0, 0, 0, 0, 1e-2]; % Variances for covariance_matrix matrix
% delta_array = [1, 0.9, 0.8];      % Normalized frequencies
% [state_transition_matrix, covariance_matrix] = get_discrete_wiener_model(L, ...
%       M, sampling_interval, sigma, delta);
%
% References:
%   - [1] Florindo, Rodrigo de Lima, Antreich, Felix, "Multi-Frequency Kalman 
%     Filter Carrier Phase Tracking for Ionospheric Scintillation Mitigation 
%     and Monitoring," Proceedings of the 37th International Technical 
%     Meeting of the Satellite Division of The Institute of Navigation 
%     (ION GNSS+ 2024), Baltimore, Maryland, September 2024, pp. 3611-3625. 
%     https://doi.org/10.33012/2024.19899
%   - [2] F. Fohlmeister, GNSS Carrier Phase Tracking under Ionospheric 
%     Scintillations, 2021-04. URL: https://elib.dlr.de/185344/.
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 ORCID: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

% Validate inputs
validateattributes(L, {'numeric'}, ...
    {'scalar', 'positive', 'integer', 'finite', 'nonnan'}, ...
    'get_discrete_wiener_model', 'L');
validateattributes(M, {'numeric'}, ...
    {'scalar', 'positive', 'integer', 'finite', 'nonnan'}, ...
    'get_discrete_wiener_model', 'M');
validateattributes(sampling_interval, {'numeric'}, ...
    {'scalar', 'positive', 'finite', 'nonnan'}, ...
    'get_discrete_wiener_model', 'sampling_interval');
validateattributes(sigma_array, {'numeric'}, ...
    {'vector', 'numel', L + M - 1, 'nonnegative', 'finite', 'nonnan'}, ...
    'get_discrete_wiener_model', 'sigma');
validateattributes(delta_array, {'numeric'}, ...
    {'vector', 'numel', L, 'positive', 'finite', 'nonnan'}, ...
    'get_discrete_wiener_model', 'delta');

% Construct F1: Transition dynamics between carriers
F1 = [zeros(L-1, 1), 2 * pi * delta_array(1:L-1).', zeros(L-1, M-2)];

% Construct F2: Dynamics for the Wiener process components
F2 = [[zeros(M-1, 1), eye(M-1)]; zeros(1, M)];
F2(1, 2) = 2 * pi * delta_array(L);

% Combine F1 and F2 into the full state transition matrix (Fw)
Fw = [[zeros(L-1, L-1), F1]; [zeros(M, L-1), F2]];

% Compute State Transition Matrix (state_transition_matrix)
state_transition_matrix = expm(Fw * sampling_interval);

% Calculate Covariance Matrix (covariance_matrix)
Q_xi = diag(sigma_array);
% Equation (15) of [1]:
% $$\int^{T_I}_{0}{e^{\mathbf{F}_w(T_I-\tau)} \mathbf{Q_{\boldsymbol{\xi}}}
% (e^{\mathbf{F}_w(T_I-\tau)})^T d\tau}$$
state_transition_exp = expm(Fw * (sampling_interval - sym('tau', 'real')));

% Compute covariance_matrix: Noise covariance matrix via symbolic integration
covariance_matrix = int(state_transition_exp * Q_xi * state_transition_exp', ...
    sym('tau', 'real'), 0, sampling_interval);

% Convert to Double
state_transition_matrix = double(state_transition_matrix);
covariance_matrix = double(covariance_matrix);
end
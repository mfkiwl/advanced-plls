function [FD, QD] = get_discrete_wiener(L, M, sampling_interval, sigma, delta)
% get_discrete_wiener
% Generates the state transition matrix and noise covariance matrix for
% line-of-sight (LOS) dynamics in a multi-frequency system.
%
% Syntax:
%   [FD, QD] = get_discrete_wiener(L, M, sampling_interval, sigma, delta)
%
% Description:
%   This function computes the state transition matrix (FD) and the noise
%   covariance matrix (QD) for a generic multi-frequency system based on
%   the mathematical framework described in Section 2.1 of the author's
%   ION GNSS+ 2024 paper.
%
% Inputs:
%   L                 - Number of carriers included in the state-space model.
%   M                 - Order of the Wiener process, indicating the degree
%                       of dynamics considered in the model.
%   sampling_interval - Time interval between samples (in seconds).
%   sigma             - Vector containing variances for the noise covariance
%                       matrix. Each variance reflects uncertainty related
%                       to clock drift and LOS dynamics for both satellite
%                       and receiver.
%   delta             - Vector of normalized frequencies for each carrier
%                       relative to the reference carrier. For example:
%                       delta = [f_L1/f_L1, f_L2/f_L1, f_L5/f_L1] for 
%                       multi-frequency tracking, or delta = 1 for 
%                       single-frequency tracking.
%
% Outputs:
%   FD                - State transition matrix of size (M+L-1) x (M+L-1),
%                       representing the dynamics of the state over time.
%   QD                - Noise covariance matrix of size (M+L-1) x (M+L-1),
%                       capturing the uncertainty in the system's state.
%
% Notes:
%   - The state transition matrix FD is derived using the matrix exponential
%     and captures the discrete-time dynamics of the system.
%   - The noise covariance matrix QD is computed through the integral of the
%     state transition matrix over time and incorporates the Wiener process.
%
% Examples:
%   % Example 1: Generate FD and QD for a 3-carrier system with a second-order
%   % Wiener process, 0.01-second sampling, variances [0, 0, 0, 0, 1e-2], and
%   % delta = [1, 0.9, 0.8].
%   L = 3; M = 3; 
%   sampling_interval = 0.01; 
%   sigma = [0, 0, 0, 0, 1e-2]; 
%   delta = [1, 0.9, 0.8];
%   [FD, QD] = get_discrete_wiener(L, M, sampling_interval, sigma, delta);
%
% References:
%   - Florindo, Rodrigo de Lima, Antreich, Felix, "Multi-Frequency Kalman 
%     Filter Carrier Phase Tracking for Ionospheric Scintillation Mitigation 
%     and Monitoring," Proceedings of the 37th International Technical 
%     Meeting of the Satellite Division of The Institute of Navigation 
%     (ION GNSS+ 2024), Baltimore, Maryland, September 2024, pp. 3611-3625. 
%     https://doi.org/10.33012/2024.19899
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 ORCID: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

% Construct State Transition Base Matrices
F1 = [zeros(L-1, 1), 2 * pi * delta(1:L-1).', zeros(L-1, M-2)];
F2 = [[zeros(M-1, 1), eye(M-1)]; zeros(1, M)];
F2(1, 2) = 2 * pi * delta(L);
Fw = [[zeros(L-1, L-1), F1]; [zeros(M, L-1), F2]];

% Compute State Transition Matrix (FD)
FD = expm(Fw * sampling_interval);

% Calculate Covariance Matrix (QD)
Q_xi = diag(sigma);
% Equation (15):
% $$\int^{T_I}_{0}{e^{\mathbf{F}_w(T_I-\tau)} \mathbf{Q_{\boldsymbol{\xi}}}
% (e^{\mathbf{F}_w(T_I-\tau)})^T d\tau}$$
state_transition_exp = expm(Fw * (sampling_interval - sym('tau', 'real')));
QD = int(state_transition_exp * Q_xi * state_transition_exp', ...
    sym('tau', 'real'), 0, sampling_interval);

% Convert to Double
FD = double(FD);
QD = double(QD);
end
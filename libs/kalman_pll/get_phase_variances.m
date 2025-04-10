function sigma2_array = get_phase_variances(C_over_N0_array_dBHz, sampling_interval)
% get_phase_variances
%
% Syntax:
%   sigma2_array = get_phase_variances(C_over_N0_array_dBHz, sampling_interval)
%
% Description:
%   Computes the phase noise variance (sigma^2) in radians squared based on the
%   carrier-to-noise density ratio (C/N0) in dB-Hz and the sampling interval in
%   seconds. This function implements the phase noise variance calculation for an
%   atan2 discriminator as given in equation 22 of [1].
%
% Inputs:
%   C_over_N0_array_dBHz - Array of carrier-to-noise density ratios (dB-Hz).
%   sampling_interval    - Sampling interval in seconds.
%
% Outputs:
%   sigma2_array         - Array of computed phase noise variances.
%
% Notes:
%   - The conversion from dB-Hz to linear scale is performed internally.
%   - This implementation is based on standard relationships between C/N0
%     and phase noise.
%   - This function can be easily used to compute the phase noise variances
%     for a multi-frequency carrier phase tracking scenario, given its
%     multivariate nature.
%
% Example:
%   % Compute phase noise variances for C/N0 values 35 and 40 dB-Hz:
%   sigma2 = get_phase_variances([35, 40], 0.01);
%
% References
% [1] R. A. M. Lopes, F. Antreich, F. Fohlmeister, M. Kriegel and H. K. Kuga,
%     "Ionospheric Scintillation Mitigation With Kalman PLLs Employing Radial
%     Basis Function Networks," in IEEE Transactions on Aerospace and
%     Electronic Systems, vol. 59, no. 5, pp. 6878-6893, Oct. 2023,
%     doi: 10.1109/TAES.2023.3281431.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Validate inputs.
    validateattributes(C_over_N0_array_dBHz, {'numeric'}, {'vector', 'real', 'finite', 'nonnan'}, mfilename, 'C_over_N0_array_dBHz');
    validateattributes(sampling_interval, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, mfilename, 'sampling_interval');
    
    % Convert dB-Hz to linear scale.
    c_over_n0_linear = 10.^(C_over_N0_array_dBHz ./ 10);
    
    % Compute variances.
    term = 1 ./ (2 * c_over_n0_linear * sampling_interval);
    sigma2_array = term .* (1 + term);
end

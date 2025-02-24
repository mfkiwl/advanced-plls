function sigma2_array = compute_phase_variances(C_over_N0_array_dBHz, sampling_interval)
    % compute_phase_variances
    % Computes the measurement noise variances based on C/N0 values.
    %
    % Syntax:
    %   sigma2_array = compute_phase_variances(C_over_N0_array_dBHz, sampling_interval)
    %
    % Description:
    %   This function calculates the phase noise variances (sigma^2) in radians squared 
    %   based on the carrier-to-noise ratio (C/N0) in dB-Hz and the sampling interval in seconds.
    %
    % Inputs:
    %   C_over_N0_array_dBHz - Array of carrier-to-noise density ratios in dB-Hz.
    %   sampling_interval    - Sampling interval in seconds.
    %
    % Outputs:
    %   sigma2_array - Array of computed measurement noise variances.
    %
    % Notes:
    %   - The calculation uses standard relationships between C/N0 and phase noise.
    %
    % Examples:
    %   sigma2_array = compute_phase_variances([35, 40], 0.01);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com
    
    C_over_N0_array_linear = 10.^(C_over_N0_array_dBHz ./ 10); % Convert dB-Hz to linear scale
    sigma2_array = (1 ./ (2 * C_over_N0_array_linear * sampling_interval)) ...
        .* (1 + 1 ./ (2 * C_over_N0_array_linear * sampling_interval)); % Compute variances
end
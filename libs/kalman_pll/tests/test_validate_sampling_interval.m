function sampling_interval = validate_sampling_interval(sampling_interval_los, sampling_interval_scint)
    % validate_sampling_interval
    % Validates consistency of sampling intervals between LOS dynamics and scintillation models.
    %
    % Syntax:
    %   sampling_interval = validate_sampling_interval(sampling_interval_los, sampling_interval_scint)
    %
    % Description:
    %   This function ensures that the sampling intervals for LOS dynamics and 
    %   scintillation models match. If they differ, an error is raised.
    %
    % Inputs:
    %   sampling_interval_los   - Sampling interval for LOS dynamics (seconds).
    %   sampling_interval_scint - Sampling interval for scintillation model (seconds).
    %
    % Outputs:
    %   sampling_interval - The validated and consistent sampling interval.
    %
    % Notes:
    %   - The function assumes both intervals are scalars and checks for equality.
    %
    % Examples:
    %   sampling_interval = validate_sampling_interval(0.01, 0.01);
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Validate consistency of sampling_interval between los and scintillation kalman_pll_config
    if sampling_interval_los ~= sampling_interval_scint
        error(['Inconsistent sampling_interval values: los (%f) vs. ' ...
               'scintillation model (%f).'], sampling_interval_los, ...
               sampling_interval_scint);
    end
    sampling_interval = sampling_interval_los; % Return consistent value
end
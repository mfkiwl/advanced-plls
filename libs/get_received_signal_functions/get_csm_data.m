function psi_csm = get_csm_data(csm_params)
% get_csm_data
%
% Syntax:
%   psi_csm = get_csm_data(csm_params)
%
% Description:
%   Generates a Cornell Scintillation Model (CSM) time series from a
%   struct of specified scintillation parameters.
%
% Inputs:
%   csm_params - Struct containing the following fields:
%       * S4 (0 <= S4 <= 1)           - Scintillation index
%       * tau0 (>0)                   - Signal intensity decorrelation time (seconds)
%       * simulation_time (>0)        - Total simulation time (seconds)
%       * sampling_interval (>0)      - Time between samples (seconds)
%
% Outputs:
%   psi_csm   - Complex field time series representing the ionospheric
%               scintillation effects.
%
% Notes:
%   - The sub-sampling parameter (Nspa = 8) can be adjusted for higher
%     accuracy (but increases computational load).
%   - The Ricean K factor is derived from the S4 index.
%   - This code uses the Cornell scintillation simulator toolkit:
%     https://gps.ece.cornell.edu/tools.php
%
% Examples:
%   csm_params = struct('S4', 0.8, 'tau0', 0.7, ...
%       'simulation_time', 600, 'sampling_interval', 0.01);
%   psi_csm = get_csm_data(csm_params);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    validateattributes(csm_params, {'struct'}, {'nonempty'}, mfilename, 'csm_params');
    
    req_fields = {'S4','tau0','simulation_time','sampling_interval'};
    for f = 1:numel(req_fields)
        if ~isfield(csm_params, req_fields{f})
            error('get_csm_data:MissingField', ...
                'csm_params must contain the field "%s".', req_fields{f});
        end
    end

    S4 = csm_params.S4;
    tau0 = csm_params.tau0;
    simulation_time = csm_params.simulation_time;
    sampling_interval = csm_params.sampling_interval;

    validateattributes(S4, {'numeric'}, ...
        {'scalar','real','finite','nonnan','>',0,'<=',1}, mfilename, 'S4');
    validateattributes(tau0, {'numeric'}, ...
        {'scalar','real','positive','finite','nonnan'}, mfilename, 'tau0');
    validateattributes(simulation_time, {'numeric'}, ...
        {'scalar','real','positive','finite','nonnan'}, mfilename, 'simulation_time');
    validateattributes(sampling_interval, {'numeric'}, ...
        {'scalar','real','positive','finite','nonnan'}, mfilename, 'sampling_interval');
    
    if simulation_time < sampling_interval
        error('get_csm_data:simulationTimeSmallerThanSamplingInterval', ...
          'The input simulation_time (%g) is smaller than the sampling_interval (%g).', ...
          simulation_time, sampling_interval);
    end

    num_samples_exact = simulation_time / sampling_interval;
    num_samples_rounded = round(num_samples_exact);

    if abs(num_samples_exact - num_samples_rounded) > eps
        warning('get_csm_data:NonIntegerRatioSamples', ...
            ['simulation_time / sampling_interval is not an integer. ' ...
             'Number of samples was rounded from %.5g to %d.'], ...
            num_samples_exact, num_samples_rounded);
    end

    Nspa = 8;

    m = max(1,1/(S4^2));
    K = sqrt(m^2 - m)/(m - sqrt(m^2 - m)); 

    [~, psi_csm, ~] = scintModel04(sampling_interval, num_samples_rounded, tau0, K, Nspa);
end

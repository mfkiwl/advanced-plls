function [psi_csm] = get_csm_data(S4,tau0,simulation_time,sampling_interval)
% get_csm_data
% Generates a Cornell Scintillation Model (CSM) time series based on
% specified scintillation parameters.
%
% Syntax:
%   [psi_csm] = get_csm_data(S4, tau0, simulation_time,sampling_interval)
%
% Description:
%   This function simulates the ionospheric scintillation complex field
%   time series based on the given scintillation parameters. The function
%   uses a modified Rician model to represent the scintillation effects
%   and calculates the required parameters such as the Ricean K factor
%   based on the specified S4 index and the decorrelation time tau0.
%
% Inputs:
%   S4              - Scintillation index (0 <= S4 <= 1), representing the
%                     severity of the amplitude fluctuations caused by scintillation.
%   tau0            - Signal intensity decorrelation time in seconds, which 
%                     describes the temporal variability of scintillation.
%   simulation_time - Total duration of the simulated time series in seconds.
%   sampling_interval 
%                   - Sampling time of the signal after integrate and dump 
%                     of a prompt correlator output
%
% Outputs:
%   psi_csm         - Complex field time series representing the
%                     ionospheric scintillation effects.
%
% Notes:
%   - The sub-sampling parameter (Nspa = 8) can be adjusted for higher 
%     accuracy at the expense of computational load.
%   - The Ricean K factor is derived based on the specified S4 index.
%
% Examples:
%   % Generate a CSM time series for S4 = 0.8, tau0 = 0.7 seconds, a
%   simulation time of 600 seconds, and a sampling time sampling_interval = 0.01:
%   [psi_csm] = get_csm_data(0.8, 0.7, 600, 0.01);
%
% References:
%   - This code uses the Cornell scintillation simulator toolkit available
%   at https://gps.ece.cornell.edu/tools.php.
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

% Input validation
validateattributes(simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_csm_data', 'simulation_time');
validateattributes(sampling_interval, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_csm_data', 'sampling_interval');
validateattributes(tau0, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_csm_data', 'tau0');
validateattributes(S4, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan', '>=', 0, '<=', 1}, 'get_csm_data', 'S4');

if simulation_time < sampling_interval
    error('get_csm_data:simulationTimeSmallerThanSamplingInterval', ...
        ['The inputed value of `simulation_time` was %g, which is smaller ' ...
        'than the value of the `sampling_interval`, %g'], simulation_time, ...
        sampling_interval)
end

% Check if simulation_time / sampling_interval is an integer
num_samples_exact = simulation_time / sampling_interval;
num_samples_rounded = round(num_samples_exact);

if abs(num_samples_exact - num_samples_rounded) > eps
    % Issue a warning if rounding was needed
    warning('get_csm_data:NonIntegerRatioSamples', ...
            ['simulation_time / sampling_interval is not an integer. ' ...
            'The number of samples was rounded from %.5g to %d.'], ...
            num_samples_exact, num_samples_rounded);
end


% Number of sub-samples that are used in calculating the averages of the 
% complex field that represents the ionospheric scintillation effect. 
% Using more sub-samples increases the accuracy of the averages
% at the expense of a greater computational burden.
Nspa = 8;

% Constant related to the scintillation intensity S4
m = max(1,1/(S4^2));

% The Ricean K (noncentrality) parameter
K = sqrt(m^2 - m)/(m - sqrt(m^2 - m));

[~,psi_csm,~] = scintModel04(sampling_interval,num_samples_rounded,tau0,K,Nspa);
end
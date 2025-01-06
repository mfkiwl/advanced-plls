function [psi_mfpsm, ps_realization] = get_mfpsm_data(S4,tau0,simulation_time,sampling_interval)
% get_mfpsm_data
% Generates Multi-Frequency Phase Screen Model (MFPSM) realizations based
% on specified scintillation parameters.
%
% Syntax:
%   [psi_mfpsm, ps_realization] = get_mfpsm_data(S4, tau0, simulation_time)
%
% Description:
%   This function calculates the phase screen realization and corresponding
%   Multi-Frequency Phase Screen Model (MFPSM) complex field time series
%   for given scintillation parameters such as the S4 index and tau0. 
%   The function also accounts for user-defined receiver location, 
%   velocity, and GPS satellite parameters.
%
% Inputs:
%   S4              - Scintillation index (0 <= S4 <= 1)
%   tau0            - Signal intensity decorrelation time in seconds
%   simulation_time - Duration of the simulation in seconds
%   sampling_interval             
%                   - Sampling time
%
% Outputs:
%   psi_mfpsm       - Output complex field timeseries (refractive + 
%   diffractive effects)
%   ps_realization  - Phase screen realization (related only to the 
%   refractive effects)
%
% Notes:
%   - This function requires manual download of broadcast files due to an 
%     unresolved issue with automatic downloads from Crustal Dynamics Data 
%     Information System (CDDIS).
%   - The receiver location defaults to Hong Kong (latitude, longitude,
%     height).
%   - Certain PRN and receiver velocity configurations may result in 
%     invalid geometry for Vdrift estimation.
%   - Original code: Main.m in the gnss-scintillation-simulator_2-param repository.
%     Modifications: Converted to a callable function for the Kalman PLL testbench.
%
% Examples:
%   % Generate MFPSM data for S4 = 0.8, tau0 = 0.7 seconds, a 
%   simulation time of 600 seconds, and a sampling time sampling_interval = 0.01:
%   [psi_mfpsm, ps_realization] = get_mfpsm_data(0.8, 0.7, 600, 0.01);
%
% References:
% - Xu D, Morton YTJ, Rino CL, Carrano CS, Jiao Y. A two-parameter 
% multifrequency GPS signal simulator for strong equatorial ionospheric 
% scintillation: modeling and parameter characterization. NAVIGATION. 2020;
% 67: 181–195. https://doi.org/10.1002/navi.350
% - GitHub fork: https://github.com/rodrigodelimaf/gnss-scintillation-simulator_2-param
% - Original repository: https://github.com/cu-sense-lab/gnss-scintillation-simulator_2-param
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com

% Input validation
validateattributes(simulation_time, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_mfpsm_data', 'simulation_time');
validateattributes(sampling_interval, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_mfpsm_data', 'sampling_interval');
validateattributes(tau0, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan'}, 'get_mfpsm_data', 'tau0');
validateattributes(S4, {'numeric'}, {'scalar', 'real', 'positive', 'finite', 'nonnan', '>=', 0, '<=', 1}, 'get_mfpsm_data', 'S4');

if simulation_time < sampling_interval
    error('get_mfpsm_data:simulationTimeSmallerThanSamplingInterval', ...
        ['The inputed value of `simulation_time` was %g, which is smaller ' ...
        'than the value of the `sampling_interval`, %g'], simulation_time, ...
        sampling_interval)
end

% Check if simulation_time / sampling_interval is an integer
num_samples_exact = simulation_time / sampling_interval;
num_samples_rounded = round(num_samples_exact);

if abs(num_samples_exact - num_samples_rounded) > eps
    % Issue a warning if rounding was needed
    warning('get_mfpsm_data:NonIntegerRatioSamples', ...
            ['simulation_time / sampling_interval is not an integer. ' ...
            'The number of samples was rounded from %.5g to %d.'], ...
            num_samples_exact, num_samples_rounded);
end

% Get the user_input settings considered the inputted S4,tau0 and
% simulation_time, as well as other fixed parameters, such as the satellite
% PRN, the receiver location and velocity, the date of the simulation, a 
% flag to generate plots regarding the simulations andthe amount of 
% frequency bands to be simulated.
user_input = get_user_input(S4,tau0,simulation_time);

% Obtain the U and rhoVeff values based on the user input S4 and tau0.
[U_mapped,rhoFVeff_mapped] = ParaMapping(user_input);

try
    satGEOM = RunPropGeomCalc(user_input,rhoFVeff_mapped);
catch ME
    error('VdriftEstimation:InvalidGeometry', ...
        ['Invalid geometry for Vdrift estimation. Try another PRN or ' ...
        'receiver velocity.\nError details: %s'], ...
        ME.message);
end

% `ps_realization` stands for phase screen realization.
% The original code was slightly modified to output the phase screen
% realization directly.
[psi_mfpsm, ps_realization] = RunGenScintFieldRealization( ...
    user_input,satGEOM,U_mapped,rhoFVeff_mapped,sampling_interval);
psi_mfpsm = psi_mfpsm(1,1:num_samples_rounded).';
ps_realization = ps_realization(1,1:num_samples_rounded).';
end

function user_input = get_user_input(S4,tau0,simulation_time)
    %%% Simulation  Settings
    %Specify the date and time in [year, moth, day, hour, minutes, seconds]
    %format
    user_input.dateTime = [2014 01 02 10 00 00];
    %Specify the simulation time with the argument of the function
    user_input.length = simulation_time;
    % Plotting figures of the simulated propagation geometry and 
    % scintillation intensity and phase? yes-1/no-0
    user_input.plotSign = 0;
    
    %%% Receiver Settings
    % Please specify receiver position as [lat(rad), lon(rad), height(m)].'
    % OBS: The actual location is at Hong Kong.
    user_input.RXPos = [0.3876 1.9942 59.6780].';
    % Please input receiver velocity as [V1,V2,V3].'
    % - V1 = east-west velocity on the earth arc (m/s, eastward +)
    % - V2 = north-south velocity on the earch arc (m/s, northward +)
    % - V3 = up-down velocity (m/s, up +)
    user_input.RXVel = [0 0 0].';
    
    %%% Satellite settings
    % Please specify satellite PRN (0~32)
    user_input.PRN = 12;
    % Please specify how many GPS frequencies to simulate
    % (1- GPS L1 only; 2 - GPS L1 and L2; 3 - GPS L1,L2, and L5)
    user_input.frequencyNo = 1;
    
    %%% Scintillation settings
    % Please specify the S4 index and tau0 for 
    % the ground observed scintillation.
    % S4 index (0~1)
    user_input.S4 = S4;
    % Signal intensity decorrelation time in sec.
    user_input.tau0 = tau0;
end
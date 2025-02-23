classdef test_get_received_signal < matlab.unittest.TestCase
% test_get_received_signal
%
% Syntax:
%   results = runtests('test_get_received_signal')
%
% Description:
%   Unit tests for the get_received_signal function. This suite tests that
%   the function simulates the baseband received signal with ionospheric
%   scintillation, thermal noise, and LOS phase dynamics correctly.
%
% Inputs:
%   (None) - No input arguments.
%
% Outputs:
%   results - Test results from MATLAB unit testing framework.
%
% Notes:
%   - Ensure that paths to the TPPSM and Cornell models are added prior to
%     testing.
%   - Intended for MATLAB R2024b.
%
% Example:
%   % Run tests from the MATLAB command window:
%   results = runtests('test_get_received_signal');
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    properties
        % Define valid parameters for testing
        C_over_N0_dBHz = 45;
        simulation_time = 300;
        settling_time = 50;
        sampling_interval = 0.01;
        S4 = 0.8;
        tau0 = 0.7;
        tppsm_scenario = 'Moderate';
        doppler_profile = [0, 1000, 0.94, 0.01];
        is_refractive_effects_removed = true;
    end
    
    methods(TestClassSetup)
        function addPaths(testCase)
            % Add the paths for the tests
            tppsm_paths = genpath(fullfile(pwd,'..','..','scintillation_models/refactored_tppsm'));
            csm_paths = genpath(fullfile(pwd,'..','..','scintillation_models/cornell_scintillation_model'));
            get_received_functions_path = fullfile(pwd,'..');
            % The paths should also include the `get_tppsm_data` function.
            all_paths = [tppsm_paths, ';', csm_paths, ';', get_received_functions_path];
            addpath(all_paths);
            testCase.addTeardown(@() rmpath(all_paths));
        end
    end
    
    methods(Test)
        function testValidNoneModel(testCase)
            % Test valid input for the 'none' scintillation model.
            [rx_sig, los_phase, psi, ps_realization] = get_received_signal(...
                testCase.C_over_N0_dBHz, 'none', testCase.doppler_profile, ...
                'simulation_time', testCase.simulation_time, ...
                'settling_time', testCase.settling_time, ...
                'sampling_interval', testCase.sampling_interval);
            
            testCase.verifyNotEmpty(rx_sig, 'Received signal should not be empty.');
            testCase.verifyNotEmpty(los_phase, 'LOS phase should not be empty.');
            % psi is expected to be a vector of ones.
            expected_size = round(testCase.simulation_time/testCase.sampling_interval);
            testCase.verifyEqual(psi, ones(expected_size, 1), 'psi must equal a ones vector for ''none'' model.');
            testCase.verifyEmpty(ps_realization, 'ps_realization should be empty for ''none'' model.');
        end
        
        function testValidCSMModel(testCase)
            % Test valid input for the 'CSM' scintillation model.
            [rx_sig, los_phase, psi, ps_realization] = get_received_signal(...
                testCase.C_over_N0_dBHz, 'CSM', testCase.doppler_profile, ...
                'S4', testCase.S4, 'tau0', testCase.tau0, ...
                'simulation_time', testCase.simulation_time, ...
                'settling_time', testCase.settling_time, ...
                'sampling_interval', testCase.sampling_interval);
            
            testCase.verifyNotEmpty(rx_sig, 'Received signal should not be empty for CSM model.');
            testCase.verifyNotEmpty(los_phase, 'LOS phase should not be empty for CSM model.');
            testCase.verifyNotEmpty(psi, 'psi should not be empty for CSM model.');
            testCase.verifyEmpty(ps_realization, 'ps_realization must be empty for CSM model.');
        end
        
        function testValidTPPSMModel(testCase)
            % Test valid input for the 'TPPSM' scintillation model.
            [rx_sig, los_phase, psi, ps_realization] = get_received_signal(...
                testCase.C_over_N0_dBHz, 'TPPSM', testCase.doppler_profile, ...
                'tppsm_scenario', testCase.tppsm_scenario, ...
                'simulation_time', testCase.simulation_time, ...
                'settling_time', testCase.settling_time, ...
                'sampling_interval', testCase.sampling_interval, ...
                'is_refractive_effects_removed', testCase.is_refractive_effects_removed);
            
            testCase.verifyNotEmpty(rx_sig, 'Received signal should not be empty for TPPSM model.');
            testCase.verifyNotEmpty(los_phase, 'LOS phase should not be empty for TPPSM model.');
            testCase.verifyNotEmpty(psi, 'psi should not be empty for TPPSM model.');
            testCase.verifyNotEmpty(ps_realization, 'ps_realization should not be empty for TPPSM model.');
        end
        
        function testMissingParametersCSM(testCase)
            % Verify error is thrown when a required parameter (tau0) for the CSM model is missing.
            testCase.verifyError(@() get_received_signal(...
                testCase.C_over_N0_dBHz, 'CSM', testCase.doppler_profile, ...
                'S4', testCase.S4, ...
                'simulation_time', testCase.simulation_time, ...
                'settling_time', testCase.settling_time), ...
                'get_received_signal:MissingParameters');
        end
        
        function testMissingParametersTPPSM(testCase)
            % Verify error is thrown when the required tppsm_scenario for the TPPSM model is missing.
            testCase.verifyError(@() get_received_signal(...
                testCase.C_over_N0_dBHz, 'TPPSM', testCase.doppler_profile, ...
                'simulation_time', testCase.simulation_time, ...
                'settling_time', testCase.settling_time), ...
                'get_received_signal:MissingParameters');
        end
        
        function testInvalidSettlingTime(testCase)
            % Verify error when the settling time exceeds the simulation time.
            testCase.verifyError(@() get_received_signal(...
                testCase.C_over_N0_dBHz, 'none', testCase.doppler_profile, ...
                'simulation_time', 100, ...
                'settling_time', 150, ...
                'sampling_interval', testCase.sampling_interval), ...
                'get_received_signal:InvalidInput');
        end
    end
end

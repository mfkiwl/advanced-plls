classdef test_get_tppsm_data < matlab.unittest.TestCase
    % test_get_tppsm_data
    %
    % Syntax:
    %   results = runtests('test_get_tppsm_data')
    %
    % Description:
    %   Unit tests for the get_tppsm_data function. This suite tests that the
    %   TPPSM simulation returns outputs with correct type and size for valid
    %   inputs, and that it errors for invalid inputs. In addition, it verifies
    %   that when an external 'rhof_veff_ratio_L1' is provided, the function uses it
    %   (i.e. by checking that the corresponding message is printed).
    %
    % Inputs:
    %   (None) - No input arguments.
    %
    % Outputs:
    %   results - Test results from MATLAB unit testing framework.
    %
    % Notes:
    %   - Ensure the correct paths to the refactored TPPSM model are added.
    %   - Intended for MATLAB R2024b.
    %
    % Example:
    %   % Run tests from the MATLAB command window:
    %   results = runtests('test_get_tppsm_data');
    %
    % Author:
    %   Rodrigo de Lima Florindo
    %   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com

    properties
        valid_scenario = 'Moderate';
        simulation_time = 300;
        sampling_interval = 0.01;
        seed = 42;
    end
    
    methods (TestClassSetup)
        function addPaths(testCase)
            % Add the paths for the tests
            tppsm_paths = genpath(fullfile(pwd,'..','..','scintillation_models','refactored_tppsm'));
            get_received_functions_path = fullfile(pwd,'..');
            all_paths = [tppsm_paths, ';' , get_received_functions_path];
            addpath(all_paths);
            % Ensure the path is removed after tests.
            testCase.addTeardown(@() rmpath(all_paths));
        end
    end
    
    methods (Test)
        function test_valid_output(testCase)
            [psi_tppsm, ps_realization, ~, ~, seedOut] = ...
                get_tppsm_data(testCase.valid_scenario, ...
                'simulation_time', testCase.simulation_time, ...
                'sampling_interval', testCase.sampling_interval, ...
                'seed', testCase.seed);
            
            % Verify outputs are non-empty and of type double.
            testCase.verifyNotEmpty(psi_tppsm);
            testCase.verifyNotEmpty(ps_realization);
            testCase.verifyClass(psi_tppsm, 'double');
            testCase.verifyClass(ps_realization, 'double');
            
            % Verify the number of samples is correct.
            expectedSamples = round(testCase.simulation_time / testCase.sampling_interval);
            testCase.verifySize(psi_tppsm, [expectedSamples, 1]);
            testCase.verifySize(ps_realization, [expectedSamples, 1]);
            
            % Verify the returned seed matches the input seed.
            testCase.verifyEqual(seedOut, testCase.seed);
        end
        
        function test_invalid_scenario(testCase)
            % Test that providing an invalid scenario produces an error.
            testCase.verifyError(@() get_tppsm_data('InvalidScenario'), ...
                'MATLAB:get_tppsm_data:unrecognizedStringChoice');
        end
        
        function test_simulation_time_less_than_sampling(testCase)
            % Test that simulation_time less than sampling_interval triggers an error.
            testCase.verifyError(@() get_tppsm_data(testCase.valid_scenario, ...
                'simulation_time', 0.005, 'sampling_interval', 0.01), ...
                'get_tppsm_data:simulationTimeSmallerThanSamplingInterval');
        end
        
        function test_non_numeric_simulation_time(testCase)
            % Test that a non-numeric simulation_time input triggers the built-in error.
            testCase.verifyError(@() get_tppsm_data(testCase.valid_scenario, 'simulation_time', 'invalid'), ...
                'MATLAB:invalidType');
        end
        
        function test_non_numeric_sampling_interval(testCase)
            % Test that a non-numeric sampling_interval input triggers the built-in error.
            testCase.verifyError(@() get_tppsm_data(testCase.valid_scenario, 'sampling_interval', 'invalid'), ...
                'MATLAB:invalidType');
        end
        
        function test_external_rhof_veff_ratio(testCase)
            % Test that providing an external rhof_veff_ratio_L1 works as expected.
            external_rhof_veff_ratio = 1;

            % Check that outputs are valid.
            [psi_tppsm, ps_realization, ~, ~, seedOut] = get_tppsm_data(testCase.valid_scenario, ...
                'simulation_time', testCase.simulation_time, ...
                'sampling_interval', testCase.sampling_interval, ...
                'seed', testCase.seed, ...
                'rhof_veff_ratio', external_rhof_veff_ratio);
            testCase.verifyNotEmpty(psi_tppsm);
            testCase.verifyNotEmpty(ps_realization);
            testCase.verifyEqual(seedOut, testCase.seed);
        end
        function test_is_enable_cmd_print_functionalty(testCase)
            % Test that providing an external rhof_veff_ratio_L1 works as expected.
            external_rhof_veff_ratio = 1;
            % Capture command-window output.
            outputStr_when_true = evalc('[psi_tppsm, ps_realization, gp, ip, seedOut] = get_tppsm_data(testCase.valid_scenario, ''is_enable_cmd_print'', true, ''simulation_time'', testCase.simulation_time, ''sampling_interval'', testCase.sampling_interval, ''seed'', testCase.seed, ''rhof_veff_ratio'', external_rhof_veff_ratio);');
            outputStr_when_false = evalc('[psi_tppsm, ps_realization, gp, ip, seedOut] = get_tppsm_data(testCase.valid_scenario, ''is_enable_cmd_print'', false, ''simulation_time'', testCase.simulation_time, ''sampling_interval'', testCase.sampling_interval, ''seed'', testCase.seed, ''rhof_veff_ratio'', external_rhof_veff_ratio);');

            % Verify that the printed message indicates the external ratio is used.
            expectedMsg = sprintf('Using externally provided rhof_veff_ratio: %g', external_rhof_veff_ratio);
            testCase.verifyNotEmpty(strfind(outputStr_when_true, expectedMsg));
            testCase.verifyEmpty(outputStr_when_false);
            testCase.verifyError(@()get_tppsm_data(testCase.valid_scenario, ...
                'is_enable_cmd_print', 'invalid'), 'MATLAB:invalidType');
        end
    end
end

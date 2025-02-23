classdef test_get_tppsm_data < matlab.unittest.TestCase
    properties
        validScenario = 'Moderate';
        simulation_time = 300;
        sampling_interval = 0.01;
        seed = 42;
    end
    
    methods (TestClassSetup)
        function addPaths(testCase)
            % Add the paths for the tests
            tppsm_paths = genpath(fullfile(pwd,'..','..','scintillation_models/refactored_tppsm'));
            get_received_functions_path = fullfile(pwd,'..');
            all_paths = [tppsm_paths, ';' , get_received_functions_path];
            addpath(all_paths);
            % Ensure the path is removed after tests.
            testCase.addTeardown(@() rmpath(all_paths));
        end
    end
    
    methods (Test)
        function testValidOutput(testCase)
            [psi_tppsm, ps_realization, general_params, irr_params_set, seedOut] = ...
                get_tppsm_data(testCase.validScenario, ...
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
            testCase.verifySize(psi_tppsm, [1, expectedSamples]);
            testCase.verifySize(ps_realization, [1, expectedSamples]);
            
            % Verify the returned seed matches the input seed.
            testCase.verifyEqual(seedOut, testCase.seed);
        end
        
        function testInvalidScenario(testCase)
            % Test that providing an invalid scenario produces an error.
            testCase.verifyError(@() get_tppsm_data('InvalidScenario'), ...
                'MATLAB:get_tppsm_data:unrecognizedStringChoice');
        end
        
        function testSimulationTimeLessThanSampling(testCase)
            % Test that simulation_time less than sampling_interval triggers an error.
            testCase.verifyError(@() get_tppsm_data(testCase.validScenario, ...
                'simulation_time', 0.005, 'sampling_interval', 0.01), ...
                'get_tppsm_data:simulationTimeSmallerThanSamplingInterval');
        end
        
        function testNonNumericSimulationTime(testCase)
            % Test that a non-numeric simulation_time input triggers the built-in error.
            testCase.verifyError(@() get_tppsm_data(testCase.validScenario, 'simulation_time', 'invalid'), ...
                'MATLAB:invalidType');
        end
        
        function testNonNumericSamplingInterval(testCase)
            % Test that a non-numeric sampling_interval input triggers the built-in error.
            testCase.verifyError(@() get_tppsm_data(testCase.validScenario, 'sampling_interval', 'invalid'), ...
                'MATLAB:invalidType');
        end
    end
end

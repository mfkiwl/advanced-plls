classdef test_get_tppsm_data < matlab.unittest.TestCase
    properties
        validScenario = 'Moderate';
        simulation_time = 300;
        sampling_interval = 0.01;
        seed = 42;
    end
    
    methods (TestClassSetup)
        function addPaths(testCase)
            % Compute repository root based on the current file location.
            % The test file is located at:
            % E:\Github\kalman_pll_testbench\libs\get_received_signal_functions\tests
            % Therefore, the repository root is four levels up.
            currentFile = mfilename('fullpath');
            repoRoot = fileparts(fileparts(fileparts(fileparts(currentFile))));
            % Define the path to the refactored TPPSM submodule.
            submodulePath = fullfile(repoRoot, 'libs', 'scintillation_models', 'refactored_tppsm');
            % Add the submodule path and all its subdirectories.
            p = genpath(submodulePath);
            addpath(p);
            % Ensure the path is removed after tests.
            testCase.addTeardown(@() rmpath(p));
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

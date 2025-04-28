classdef test_build_extended_kf < matlab.unittest.TestCase
    properties
        F_los
        Q_los
        general_config
        sampling_interval
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            % Create a simple LOS dynamics example.
            testCase.F_los = eye(2);
            testCase.Q_los = eye(2);
            testCase.sampling_interval = 0.01;
            % Set C_over_N0_array_dBHz. For testing, it's a scalar.
            testCase.general_config.C_over_N0_array_dBHz = 35;
        end
    end
    
    methods(Test)
        function testNoAugmentation(testCase)
            % --- Test the 'none' branch for build_extended_kf ---
            aug_data.augmentation_type = 'none';
            
            % Call the function under test.
            [F, Q, H_handlej, R, W] = build_extended_kf(...
                testCase.F_los, testCase.Q_los, aug_data, ...
                testCase.general_config, testCase.sampling_interval);
            
            % For the 'none' branch, F and Q should match the LOS matrices.
            testCase.verifyEqual(F, testCase.F_los, ...
                'F should equal F_los when no augmentation is applied.');
            testCase.verifyEqual(Q, testCase.Q_los, ...
                'Q should equal Q_los when no augmentation is applied.');
            testCase.verifyEqual(W, zeros(size(testCase.F_los,1), 1), ...
                'W should be a zero vector with length equal to the number of KF states.');
            
            % Compute the expected R dynamically from get_phase_variances.
            % Note: This makes the test independent of an assumption on the value.
            variances = get_phase_variances(testCase.general_config.C_over_N0_array_dBHz, testCase.sampling_interval);
            expectedR = diag(repelem(variances,2) / 2);
            testCase.verifyEqual(R, expectedR, ...
                'R must match the computed I_Q_separate_variances from get_phase_variances.');
            
            % Test the measurement Jacobian handle.
            % For the 'none' branch the measurement model is:
            %   h(x) = [ external_amplitude*cos(phase);
            %            external_amplitude*sin(phase) ]
            % with the KF state vector x containing the phase in its first element.
            external_amplitude = 10;
            phase = pi/6;  % 30 degrees in radians.
            % Create a KF state vector (the length is determined by F_los; here length 2).
            x = [phase; 0]; 
            
            % Call the function handle.
            Hj = H_handlej(external_amplitude, x);
            
            % Expected derivatives:
            %   dI/dphase = -external_amplitude*sin(phase)
            %   dQ/dphase =  external_amplitude*cos(phase)
            expected_dI_dphase = -external_amplitude * sin(phase);
            expected_dQ_dphase = external_amplitude * cos(phase);
            expectedHj = zeros(2, length(x));
            expectedHj(:, 1) = [expected_dI_dphase; expected_dQ_dphase];
            
            testCase.verifyEqual(Hj, expectedHj, 'AbsTol', 1e-4, ...
                'The measurement Jacobian does not match the expected values.');
        end
        
        function testUnsupportedAugmentation(testCase)
            % --- Test that an unsupported augmentation type produces an error ---
            aug_data.augmentation_type = 'arfit';
            
            % Define an anonymous function for calling build_extended_kf.
            testFunc = @() build_extended_kf(...
                testCase.F_los, testCase.Q_los, aug_data, ...
                testCase.general_config, testCase.sampling_interval);
            
            testCase.verifyError(testFunc, 'MATLAB:UndefinedFunction', ...
                'An unsupported augmentation type should trigger the proper error.');
        end
    end
end

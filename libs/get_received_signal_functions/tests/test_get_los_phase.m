classdef test_get_los_phase < matlab.unittest.TestCase
% test_get_los_phase
% Unit tests for the `get_los_phase` function, which computes the line-of-sight 
% (LOS) phase time series based on a Taylor series expansion defined by the
% `doppler_profile` input.
%
% This test suite validates the correctness and robustness of the `get_los_phase`
% function by testing its behavior under various scenarios, including:
%
% **Functional Tests**:
%   - Validate the correct progression of the LOS phase time series.
%   - Ensure numerical stability in the computation of LOS phase.
%
% **Edge Case Tests**:
%   1. `simulation_time == sampling_interval`: No warnings or errors expected.
%   2. `simulation_time < sampling_interval`: Function must throw an error.
%   3. `simulation_time == 0`: Function must throw an error.
%   4. `sampling_interval == 0`: Function must throw an error.
%   5. Non-integer `simulation_time / sampling_interval` ratio:
%       - Validate that a warning is issued.
%       - Ensure rounding behavior is correct for such cases.
%
% **Validation Tests**:
%   - Dynamically validate error handling for invalid inputs, including:
%       1. Negative or zero values for `simulation_time`, `sampling_interval`,
%          or invalid `doppler_profile` values.
%       2. Non-numeric, complex, or non-scalar values for the parameters.
%       3. Boundary conditions, such as empty or infinite inputs.
%
% Example:
%   To run the test suite:
%       results = runtests('test_get_los_phase');
%       disp(results);
%
%   To debug specific test cases:
%       results = runtests('test_get_los_phase', 'ProcedureName', 'test_invalid_inputs');
%       disp(results);
%
% Author 1: Rodrigo de Lima Florindo
% Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
% Author's 1 Email: rdlfresearch@gmail.com
properties
    % Test parameters (valid inputs)
    simulation_time = 300; % seconds
    sampling_interval = 0.01; % seconds
    doppler_profile = [0, 1000, 0.94]; % Initial phase, Doppler shift, Doppler drift

    % Invalid input cases, initialized in the TestClassSetup
    invalid_inputs
end

properties (TestParameter)
    % Parameterized property for input names
    input_name = struct(...
        'simulation_time', 'simulation_time', ...
        'sampling_interval', 'sampling_interval', ...
        'doppler_profile', 'doppler_profile' ...
    );
end

methods (TestClassSetup)
    function classSetup(testCase)
        % Add parent directory to path
        pathToAdd = fullfile(pwd, '..');
        if ~contains(path, pathToAdd)
            addpath(pathToAdd);
            testCase.addTeardown(@() rmpath(pathToAdd));
        end

        % Define invalid input cases as cell arrays to avoid dimension mismatch
        testCase.invalid_inputs = struct(...
            'simulation_time', { ...
                { 'invalid', 'MATLAB:get_los_phase:invalidType'; ...
                  true, 'MATLAB:get_los_phase:invalidType'; ...
                  {300}, 'MATLAB:get_los_phase:invalidType'; ...
                  [300, 400], 'MATLAB:get_los_phase:expectedScalar'; ...
                  zeros(3, 3), 'MATLAB:get_los_phase:expectedScalar'; ...
                  -300, 'MATLAB:get_los_phase:expectedPositive'; ...
                  0, 'MATLAB:get_los_phase:expectedPositive'; ...
                  Inf, 'MATLAB:get_los_phase:expectedFinite'; ...
                  -Inf, 'MATLAB:get_los_phase:expectedPositive'; ...
                  NaN, 'MATLAB:get_los_phase:expectedFinite'; ...
                  [], 'MATLAB:get_los_phase:expectedScalar'; ...
                  300 + 1j, 'MATLAB:get_los_phase:expectedReal' ...
                } ...
            }, ...
            'sampling_interval', { ...
                { 'invalid', 'MATLAB:get_los_phase:invalidType'; ...
                  true, 'MATLAB:get_los_phase:invalidType'; ...
                  {0.01}, 'MATLAB:get_los_phase:invalidType'; ...
                  [0.01, 0.02], 'MATLAB:get_los_phase:expectedScalar'; ...
                  zeros(3, 3), 'MATLAB:get_los_phase:expectedScalar'; ...
                  -0.01, 'MATLAB:get_los_phase:expectedPositive'; ...
                  0, 'MATLAB:get_los_phase:expectedPositive'; ...
                  Inf, 'MATLAB:get_los_phase:expectedFinite'; ...
                  -Inf, 'MATLAB:get_los_phase:expectedPositive'; ...
                  NaN, 'MATLAB:get_los_phase:expectedFinite'; ...
                  [], 'MATLAB:get_los_phase:expectedScalar'; ...
                  0.01 + 1j, 'MATLAB:get_los_phase:expectedReal' ...
                } ...
            }, ...
            'doppler_profile', { ...
                { NaN, 'MATLAB:get_los_phase:expectedFinite'; ...
                  Inf, 'MATLAB:get_los_phase:expectedFinite'; ...
                  -Inf, 'MATLAB:get_los_phase:expectedFinite'; ...
                  [0, NaN], 'MATLAB:get_los_phase:expectedFinite'; ...
                  [0, Inf], 'MATLAB:get_los_phase:expectedFinite'; ...
                  [0, -Inf], 'MATLAB:get_los_phase:expectedFinite'; ...
                  [], 'MATLAB:get_los_phase:expectedRow'; ...
                  magic(3), 'MATLAB:get_los_phase:expectedRow'; ...
                  [0;1000;0.94], 'MATLAB:get_los_phase:expectedRow' ...
                } ...
            } ...
        );
    end
end

methods (Test)
    %% Edge Case Scenarios

    %%% Test Case: simulation_time == sampling_interval
    function test_simulation_time_equals_sampling_interval(testCase)
        % simulation_time = sampling_interval; expect no warning
        los_phase = get_los_phase(0.01, 0.01, testCase.doppler_profile);
        testCase.verifySize(los_phase, [1, 1], ...
            'LOS phase output does not have the expected size for a single sample.');
    end

    %%% Test Case: simulation_time < sampling_interval
    function test_simulation_time_smaller_than_sampling_interval(testCase)
        % simulation_time < sampling_interval; expect error
        testCase.verifyError(@() get_los_phase(0.001, 0.01, testCase.doppler_profile), ...
                             'get_los_phase:simulationTimeSmallerThanSamplingInterval');
    end

    %%% Test Case: simulation_time == 0
    function test_simulation_time_zero(testCase)
        % simulation_time = 0; expect error
        testCase.verifyError(@() get_los_phase(0, testCase.sampling_interval, testCase.doppler_profile), ...
                             'MATLAB:get_los_phase:expectedPositive');
    end

    %%% Test Case: sampling_interval == 0
    function test_sampling_interval_zero(testCase)
        % sampling_interval = 0; expect error
        testCase.verifyError(@() get_los_phase(testCase.simulation_time, 0, testCase.doppler_profile), ...
                             'MATLAB:get_los_phase:expectedPositive');
    end

    %%% Test Case: simulation_time / sampling_interval non-integer
    function test_simulation_time_div_sampling_interval_non_integer(testCase)
        % simulation_time / sampling_interval is non-integer; expect a warning
        testCase.verifyWarning(@() get_los_phase(0.011, 0.01, testCase.doppler_profile), ...
                               'get_los_phase:NonIntegerRatio');
    end
    %%% Test Case: Zero Doppler Profile
    function test_zero_doppler_profile(testCase)
        % Test with doppler_profile = 0
        los_phase = get_los_phase(testCase.simulation_time, testCase.sampling_interval, 0);
        testCase.verifyEqual(los_phase, zeros(size(los_phase)), ...
            'LOS phase should be all zeros for doppler_profile = [0].');
    end

    %%% Test Case: High-Order Doppler Profile
    function test_high_order_doppler_profile(testCase)
        % Test with higher-order coefficients in doppler_profile
        los_phase = get_los_phase(testCase.simulation_time, testCase.sampling_interval, [0, 1000, 0.94, 0.1, 0.001]);
        num_samples = round(testCase.simulation_time / testCase.sampling_interval);
        testCase.verifySize(los_phase, [num_samples, 1], ...
            'LOS phase output does not have the expected size for a high-order doppler_profile.');
    end
    %% Parameterized Validation Tests
    function test_invalid_inputs(testCase, input_name)
        testCase.run_validation_tests(input_name, testCase.invalid_inputs.(input_name));
    end
end

methods
    function run_validation_tests(testCase, input_name, invalidCases)
        % Iterate through invalid cases and verify errors
        for invalidCase = invalidCases.'
            % Extract invalid input and expected error
            invalidInput = invalidCase{1, 1};
            expectedError = invalidCase{2, 1};

            % Generate function inputs and convert invalid input to string for messages
            functionInputs = testCase.generate_inputs(input_name, invalidInput);
            invalidInputStr = testCase.safe_input_strings(invalidInput);

            % Run the validation test and check for expected error
            testCase.verifyError(@() get_los_phase(functionInputs{:}), expectedError, ...
                sprintf('Validation failed for %s with input: %s', input_name, invalidInputStr));
        end
    end

    function inputs = generate_inputs(testCase, fieldName, value)
        % Create inputs for the get_los_phase function
        inputs = {...
            testCase.simulation_time, ...
            testCase.sampling_interval, ...
            testCase.doppler_profile};
        if strcmp(fieldName, 'doppler_profile')
            inputs{3} = value;
        elseif strcmp(fieldName, 'simulation_time')
            inputs{1} = value;
        elseif strcmp(fieldName, 'sampling_interval')
            inputs{2} = value;
        end
    end

    function str = safe_input_strings(~, input)
        % Convert input to string safely for error messages
        if ischar(input) || isstring(input)
            str = char(input);
        elseif isnumeric(input) || islogical(input)
            str = mat2str(input);
        else
            str = '<unconvertible input>';
        end
    end
end
end

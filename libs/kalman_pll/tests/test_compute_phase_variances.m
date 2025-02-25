classdef test_compute_phase_variances < matlab.unittest.TestCase
% test_compute_phase_variances
%
% Syntax:
%   results = runtests('test_compute_phase_variances')
%
% Description:
%   Unit tests for the compute_phase_variances function. This test suite
%   verifies that the function computes the expected phase noise variances
%   based on provided C/N0 values (in dB-Hz) and a sampling interval (in seconds).
%   It tests the output size, positivity, approximate numerical values, and
%   that invalid inputs trigger appropriate errors.
%
% Example:
%   % Run the test suite:
%   results = runtests('test_compute_phase_variances');
%   disp(results);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    methods(TestClassSetup)
        function add_parent_path(test_case)
            % add_parent_path - Add the parent directory containing compute_phase_variances.
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            addpath(parent_dir);
            test_case.addTeardown(@() rmpath(parent_dir));
        end
    end

    methods(Test)
        function test_output_values(test_case)
            % Test that computed phase noise variances are within tolerance.
            c_over_n0_array_dbhz = [35, 40];
            sampling_interval = 0.01;
            sigma2_array = compute_phase_variances(c_over_n0_array_dbhz, sampling_interval);
            
            % Expected approximate values:
            % For 35 dB-Hz: linear ~ 3162.28 => term ~ 1/(2*3162.28*0.01)=0.01582, sigma2 ~ 0.01582*(1+0.01582) ≈ 0.01607
            % For 40 dB-Hz: linear ~ 10000   => term = 1/(2*10000*0.01)=0.005, sigma2 ~ 0.005*(1+0.005) ≈ 0.005025
            expected_values = [0.01607, 0.005025];
            tol = 1e-4;
            test_case.verifySize(sigma2_array, size(c_over_n0_array_dbhz), 'Output size mismatch.');
            test_case.verifyLessThanOrEqual(abs(sigma2_array - expected_values), tol, ...
                sprintf('Computed sigma2 values are not within tolerance %.5g.', tol));
        end

        function test_output_positivity(test_case)
            % Verify that the output variances are positive.
            c_over_n0_array_dbhz = [30, 35, 40];
            sampling_interval = 0.01;
            sigma2_array = compute_phase_variances(c_over_n0_array_dbhz, sampling_interval);
            test_case.verifyGreaterThan(sigma2_array, 0, 'All computed variances should be positive.');
        end

        function test_invalid_sampling_interval_error(test_case)
            % Test that a non-numeric sampling_interval triggers an error.
            c_over_n0_array_dbhz = [35, 40];
            invalid_sampling_interval = 'invalid';
            test_case.verifyError(@() compute_phase_variances(c_over_n0_array_dbhz, invalid_sampling_interval), ...
                'MATLAB:compute_phase_variances:invalidType');
        end

        function test_invalid_c_over_n0_error(test_case)
            % Test that a non-numeric C_over_N0_array_dBHz triggers an error.
            invalid_c_over_n0 = 'invalid';
            sampling_interval = 0.01;
            test_case.verifyError(@() compute_phase_variances(invalid_c_over_n0, sampling_interval), ...
                'MATLAB:compute_phase_variances:invalidType');
        end
    end
end

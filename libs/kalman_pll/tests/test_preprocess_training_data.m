classdef test_preprocess_training_data < matlab.unittest.TestCase
% test_preprocess_training_data
%
% Syntax:
%   results = runtests('test_preprocess_training_data')
%
% Description:
%   Unit tests for the preprocess_training_data function, which
%   generates phase data from either CSM, TPPSM, or "none" scintillation
%   models, optionally removing refractive effects for TPPSM.
%
% Example:
%   % Run the test suite:
%   results = runtests('test_preprocess_training_data');
%   disp(results);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
    %   Email: rdlfresearch@gmail.com

    methods(TestClassSetup)
        function add_parent_path(test_case)
            % add_parent_path - Add the parent directory of the test folder 
            % so that preprocess_training_data.m and dependencies can be found.
            parent_dir = fileparts(fileparts(mfilename('fullpath')));
            get_received_signal_functions_dir = [fileparts(parent_dir),'\get_received_signal_functions'];
            tppsm_paths = genpath([fileparts(parent_dir),'\scintillation_models\refactored_tppsm']);
            csm_paths = genpath([fileparts(parent_dir),'\scintillation_models\cornell_scintillation_model']);
            addpath(parent_dir);
            addpath(get_received_signal_functions_dir);
            addpath(tppsm_paths);
            addpath(csm_paths);

            test_case.addTeardown(@() rmpath(parent_dir,get_received_signal_functions_dir, ...
                tppsm_paths, csm_paths));
        end
    end

    methods(Test)
        function test_valid_csm_config(test_case)
            % test_valid_csm_config
            % Ensures a valid CSM config produces nonempty training_data 
            % without errors.
            
            csm_config = struct( ...
                'scintillation_model', 'CSM', ...
                'S4', 0.8, ...
                'tau0', 0.7, ...
                'simulation_time', 600, ...
                'sampling_interval', 0.01 ...
            );
            
            training_data = preprocess_training_data(csm_config);
            test_case.verifyNotEmpty(training_data, ...
                'Expected nonempty training_data for valid CSM config.');
        end
        
        function test_valid_tppsm_config(test_case)
            % test_valid_tppsm_config
            % Ensures a valid TPPSM config runs without error, 
            % default is_refractive_effects_removed = false if missing.
            
            tppsm_config = struct( ...
                'scintillation_model', 'TPPSM', ...
                'scenario', 'Moderate', ...
                'simulation_time', 300, ...
                'sampling_interval', 0.01 ...
            );
            
            training_data = preprocess_training_data(tppsm_config);
            test_case.verifyNotEmpty(training_data, ...
                'Expected nonempty training_data for valid TPPSM config (no refractive removal).');
        end
        
        function test_tppsm_refractive_effects_removed(test_case)
            % test_tppsm_refractive_effects_removed
            % Check if is_refractive_effects_removed = true triggers no error
            tppsm_config = struct( ...
                'scintillation_model', 'TPPSM', ...
                'scenario', 'Severe', ...
                'simulation_time', 200, ...
                'sampling_interval', 0.01, ...
                'is_refractive_effects_removed', true ...
            );
            
            training_data = preprocess_training_data(tppsm_config);
            test_case.verifyNotEmpty(training_data, ...
                'Expected nonempty training_data with refractive effects removal for TPPSM.');
        end
        
        function test_model_none(test_case)
            % test_model_none
            % If scintillation_model is 'none', training_data should be empty
            none_config = struct('scintillation_model','none');
            td = preprocess_training_data(none_config);
            test_case.verifyEmpty(td, ...
                'Expected empty training_data when scintillation_model = ''none''.');
        end
        
        function test_missing_scintillation_model(test_case)
            % test_missing_scintillation_model
            % Omitting scintillation_model field should raise an error
            invalid_config = struct('S4', 0.8, 'tau0', 0.7);
            test_case.verifyError(@() preprocess_training_data(invalid_config), ...
                'preprocess_training_data:MissingScintillationModelField', ...
                'Expected an error for missing scintillation_model field.');
        end
        
        function test_invalid_model_type(test_case)
            % test_invalid_model_type
            % scintillation_model not a char or string should raise an error
            invalid_config = struct('scintillation_model', 123);
            test_case.verifyError(@() preprocess_training_data(invalid_config), ...
                'preprocess_training_data:InvalidModelType', ...
                'Expected an error for invalid scintillation_model type.');
        end
        
        function test_missing_csm_field(test_case)
            % test_missing_csm_field
            % For CSM, if S4, tau0, simulation_time, sampling_interval are missing
            % we expect an error.
            csm_config = struct( ...
                'scintillation_model', 'CSM', ...
                'tau0', 0.7, ...
                'simulation_time', 600, ...
                'sampling_interval', 0.01 ...
            );
            % Missing S4
            test_case.verifyError(@() preprocess_training_data(csm_config), ...
                'preprocess_training_data:MissingCSMField', ...
                'Expected an error due to missing S4 in CSM config.');
        end
        
        function test_missing_tppsm_field(test_case)
            % test_missing_tppsm_field
            % For TPPSM, scenario, simulation_time, sampling_interval are required.
            tp_config = struct( ...
                'scintillation_model', 'TPPSM', ...
                'simulation_time', 300, ...
                'sampling_interval', 0.01, ...
                'is_refractive_effects_removed', false ...
            );
            % Missing 'scenario'
            test_case.verifyError(@() preprocess_training_data(tp_config), ...
                'preprocess_training_data:MissingTPPSMField', ...
                'Expected an error for missing scenario in TPPSM config.');
        end
        
        function test_struct_is_empty(test_case)
            % test_struct_is_empty
            % Passing an empty struct raises error
            empty_struct = struct([]);
            test_case.verifyError(@() preprocess_training_data(empty_struct), ...
                'MATLAB:preprocess_training_data:expectedNonempty', ...
                'Expected a mustBeNonempty error when config is an empty struct.');
        end
        
        function test_invalid_input_type(test_case)
            % test_invalid_input_type
            % Passing a non-struct input should raise an error
            invalid_input = 42;
            test_case.verifyError(@() preprocess_training_data(invalid_input), ...
                'MATLAB:preprocess_training_data:invalidType', ...
                'Expected an error for invalid (non-struct) input type.');
        end
        
        function test_unsupported_model(test_case)
            % test_unsupported_model
            % If scintillation_model is an unrecognized string, we expect an error
            invalid_config = struct('scintillation_model', 'InvalidModelXYZ');
            test_case.verifyError(@() preprocess_training_data(invalid_config), ...
                'preprocess_training_data:UnsupportedModel', ...
                'Expected an error for unsupported scintillation model.');
        end
    end
end

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
            get_received_signal_functions_dir = fullfile(fileparts(parent_dir),'get_received_signal_functions');
            tppsm_paths = genpath(fullfile(fileparts(parent_dir),'scintillation_models','refactored_tppsm'));
            csm_paths = genpath(fullfile(fileparts(parent_dir),'scintillation_models','cornell_scintillation_model'));
            addpath(parent_dir);
            addpath(genpath(get_received_signal_functions_dir));
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
                'sampling_interval', 0.01, ...
                'is_unwrapping_used', false ...
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
                'sampling_interval', 0.01, ...
                'is_refractive_effects_removed', true, ...
                'is_unwrapping_used', false ...
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
                'is_refractive_effects_removed', true, ...
                'is_unwrapping_used', false ...
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
    end
end

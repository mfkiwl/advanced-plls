function training_data = preprocess_training_data(scint_training_data_cfg)
% preprocess_training_data
%
% Syntax:
%   training_data = preprocess_training_data(scint_training_data_cfg)
%
% Description:
%   This function generates training data for a given scintillation model.
%   It extracts the phase information from the scintillation field, optionally
%   removing refractive effects for the TPPSM model.
%
% Inputs:
%   scint_training_data_cfg - Struct containing scintillation model settings.
%       For CSM, expected fields:
%           scintillation_model - Must be 'CSM'
%           S4                  - Scintillation index (0 <= S4 <= 1)
%           tau0                - Signal decorrelation time (positive scalar)
%           simulation_time     - Duration of simulation (positive scalar)
%           sampling_interval   - Sampling interval (positive scalar)
%           is_unwrapping_used - Boolean flag (true or false) to use
%                                unwrapping function for training the 
%                                augmentation models.
%
%       For TPPSM, expected fields:
%           scintillation_model - Must be 'TPPSM'
%           scenario            - A string specifying the scenario ('Weak', 'Moderate', 'Severe')
%           simulation_time     - Duration of simulation (positive scalar)
%           sampling_interval   - Sampling interval (positive scalar)
%           is_refractive_effects_removed - Boolean flag (true or false)
%           is_unwrapping_used - Boolean flag (true or false) to use
%                                unwrapping function for training the 
%                                augmentation models.
%
% Outputs:
%   training_data - Phase data extracted from the scintillation field.
%
% Notes:
%   - Calls external functions (get_csm_data or get_tppsm_data) to generate scintillation fields.
%   - For TPPSM, if is_refractive_effects_removed is true, refractive effects are removed by
%     multiplying the field with the complex conjugate of the phase screen realization.
%
% Examples:
%   % For CSM:
%   csmConfig = struct('scintillation_model', 'CSM', 'S4', 0.8, 'tau0', 0.7, ...
%                      'simulation_time', 600, 'sampling_interval', 0.01);
%   training_data = preprocess_training_data(csmConfig);
%
%   % For TPPSM:
%   tpConfig = struct('scintillation_model', 'TPPSM', 'scenario', 'Moderate', ...
%                      'simulation_time', 300, 'sampling_interval', 0.01, ...
%                      'is_refractive_effects_removed', true);
%   training_data = preprocess_training_data(tpConfig);
%
% Dependencies:
%   get_csm_data, get_tppsm_data, struct_to_nv_pairs.
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    switch upper(scint_training_data_cfg.scintillation_model)
        case 'CSM'
            % Generate CSM data
            csm_config = rmfield(scint_training_data_cfg,'scintillation_model');
            scint_complex_field = get_csm_data(csm_config);
            training_data = angle(scint_complex_field);
            if scint_training_data_cfg.is_unwrapping_used
                training_data = unwrap(training_data);
            end

        case 'TPPSM'
            % Build name-value pairs from the struct except the "scintillation_model" field
            tppsm_config = rmfield(scint_training_data_cfg, {'scintillation_model','is_refractive_effects_removed','is_unwrapping_used'});
            nv = struct_to_nv_pairs(tppsm_config);
            [scint_complex_field, ps_realization, ~, ~, ~] = ...
                get_tppsm_data(scint_training_data_cfg.scenario, nv{:});

            % Optionally remove refractive effects
            if scint_training_data_cfg.is_refractive_effects_removed
                scint_complex_field = scint_complex_field .* exp(-1j * ps_realization);
            end
            training_data = angle(scint_complex_field);
            if scint_training_data_cfg.is_unwrapping_used
                training_data = unwrap(training_data);
            end
        case 'NONE'
            % For 'none', we can define training_data as an empty array
            training_data = [];

        otherwise
            error('preprocess_training_data:UnsupportedModel', ...
                'Unsupported scintillation model: %s', model);
    end
end

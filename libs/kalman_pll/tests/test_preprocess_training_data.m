function training_data = preprocess_training_data(scintillation_training_data_config)
% preprocess_training_data
%
% Syntax:
%   training_data = preprocess_training_data(scintillation_training_data_config)
%
% Description:
%   This function generates training data for a given scintillation model.
%   It extracts the phase information from the scintillation field, optionally
%   removing refractive effects for the TPPSM model.
%
% Inputs:
%   scintillation_training_data_config - Struct containing scintillation model settings.
%       For CSM, expected fields:
%           S4              - Scintillation index (0 <= S4 <= 1)
%           tau0            - Signal decorrelation time (positive scalar)
%           simulation_time - Duration of simulation (positive scalar)
%           sampling_interval - Sampling interval (positive scalar)
%
%       For TPPSM, expected fields:
%           scintillation_model - Must be 'TPPSM'
%           scenario            - A string specifying the scenario ('Weak', 'Moderate', 'Severe')
%           simulation_time     - Duration of simulation (positive scalar)
%           sampling_interval   - Sampling interval (positive scalar)
%           is_refractive_effects_removed - Boolean flag (true or false)
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
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    model = scintillation_training_data_config.scintillation_model;
    
    switch upper(model)
        case 'CSM'
            % For CSM, extract required fields.
            csm_config = scintillation_training_data_config;
            scint_complex_field = get_csm_data(csm_config.S4, csm_config.tau0, ...
                csm_config.simulation_time, csm_config.sampling_interval);
            training_data = angle(scint_complex_field); % Extract phase data
            
        case 'TPPSM'
            % For TPPSM, extract required fields.
            tppsm_config = scintillation_training_data_config;
            % Automatically build name-value pairs from the struct (except for 'scintillation_model')
            nv = struct_to_nv_pairs(rmfield(tppsm_config, 'scintillation_model'));
            [scint_complex_field, ps_realization, ~, ~, ~] = get_tppsm_data(tppsm_config.scenario, nv{:});
            if tppsm_config.is_refractive_effects_removed
                % Remove refractive effects.
                scint_complex_field = scint_complex_field .* exp(-1j * ps_realization);
            end
            training_data = angle(scint_complex_field); % Extract phase data
            
        case 'NONE'
            training_data = [];

        otherwise
            error('preprocess_training_data:UnsupportedModel', 'Unsupported scintillation model: %s', model);
    end
end
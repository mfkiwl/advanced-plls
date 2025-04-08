function [kalman_pll_config, is_cache_used] = get_cached_kalman_pll_config(general_config, cache_file, is_enable_cmd_print)
% get_cached_kalman_pll_config 
% 
% Loads a cached Kalman PLL configuration if available and valid.
%
% Syntax:
%   [kalman_pll_config, is_cache_used] = get_cached_kalman_pll_config(general_config, cache_file, is_enable_cmd_print)
%
% Description:
%   Loads a cached Kalman PLL configuration from the specified cache file if it exists.
%   The cache is expected to contain an outer field corresponding to the scintillation model
%   and, within that, an inner structure with fields for different KF types
%   ('standard', 'extended', 'unscented', 'cubature'). If the cache is missing or the required
%   settings are absent, a new configuration is initialized.
%
% Inputs:
%   general_config - Struct containing user-defined PLL configuration settings. It must include:
%       .scintillation_training_data_config.scintillation_model
%       .kf_type (e.g., 'standard', 'extended', 'unscented', or 'cubature')
%       .is_use_cached_settings (boolean)
%
%   cache_file - String specifying the path to the cache file.
%
%   is_enable_cmd_print - Boolean flag for enabling command-line prints.
%
% Outputs:
%   kalman_pll_config - A struct with the computed or loaded Kalman filter settings. Organized as:
%                       kalman_pll_config.(scint_model).(kf_type)
%
%   is_cache_used     - Boolean flag that is true if the cached configuration was used.
%
% Example:
%   general_config = struct(...);  % User-defined configuration settings
%   [pll_cfg, cache_used] = get_cached_kalman_pll_config(general_config, 'cache/kalman_pll_cache.mat', true);
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Initialize the output flag.
    is_cache_used = false;
    
    % Extract scintillation model and KF type from the general configuration.
    scint_model = general_config.scintillation_training_data_config.scintillation_model;
    kf_type = lower(string(general_config.kf_type)); % Ensure consistent (lower-case) naming

    % Load from cache if available.
    if isfile(cache_file)
        if is_enable_cmd_print
            fprintf('Cache file found. Loading cached kalman_pll_config.\n');
        end
        load(cache_file, 'kalman_pll_config');
        
        % Check if the outer field for the scintillation model exists.
        if isfield(kalman_pll_config, scint_model)
            inner_config = kalman_pll_config.(scint_model);
            % Check if the inner configuration for the current KF type exists and is nonempty.
            if isfield(inner_config, kf_type) && ~isempty(inner_config.(kf_type))
                if general_config.is_use_cached_settings
                    if is_enable_cmd_print
                        fprintf('Using cached values for scintillation model "%s" and KF type "%s".\n', scint_model, kf_type);
                    end
                    is_cache_used = true;
                    return;
                else
                    if is_enable_cmd_print
                        fprintf('Recomputing values for scintillation model "%s" and KF type "%s" and updating cache.\n', scint_model, kf_type);
                    end
                end
            else
                if is_enable_cmd_print
                    fprintf('Cache found for "%s" but missing KF type "%s" settings. Computing and caching new values.\n', scint_model, kf_type);
                end
            end
        else
            if is_enable_cmd_print
                fprintf('Cache found but missing scintillation model "%s". Computing and caching new values.\n', scint_model);
            end
        end
    else
        if is_enable_cmd_print
            fprintf('No cache file found. Initializing kalman_pll_config.\n');
        end
        % Initialize a new kalman_pll_config structure with all expected fields.
        kalman_pll_config = struct(...
            'CSM',   struct('standard', [], 'extended', [], 'unscented', [], 'cubature', []), ...
            'TPPSM', struct('standard', [], 'extended', [], 'unscented', [], 'cubature', []), ...
            'none',    struct('standard', [], 'extended', [], 'unscented', [], 'cubature', []));
    end
end

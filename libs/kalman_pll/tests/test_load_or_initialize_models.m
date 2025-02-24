function [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file)
    % load_or_initialize_models
    % Loads cached Kalman filter configuration or initializes a new one.
    %
    % Syntax:
    %   [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file)
    %
    % Description:
    %   This function attempts to load cached Kalman filter configurations 
    %   from a specified file. If the cache file does not exist or is invalid, 
    %   it initializes a new configuration.
    %
    % Inputs:
    %   config     - Struct containing user configuration.
    %   cache_file - Path to the cache file for storing Kalman filter configurations.
    %
    % Outputs:
    %   kalman_pll_config - The Kalman filter configuration struct.
    %   is_cache_used     - Boolean flag indicating whether cached values were used.
    %
    % Notes:
    %   - If the cache file exists but does not contain the required configuration 
    %     for the specified scintillation model, a new configuration is computed.
    %
    % Examples:
    %   [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, 'kalman_pll_cache.mat');
    %
    % Author 1: Rodrigo de Lima Florindo
    % Author's 1 Orcid: https://orcid.org/0000-0003-0412-5583
    % Author's 1 Email: rdlfresearch@gmail.com

    % Initialize the output flag
    is_cache_used = false;

    % Check if the cache file exists
    if isfile(cache_file)
        % Load existing cached kalman_pll_config
        fprintf('Cache file found. Loading cached kalman_pll_config.\n');
        load(cache_file, 'kalman_pll_config');

        % Check if the specific scintillation model in the cache is non-empty
        if isfield(kalman_pll_config, config.training_scint_model) && ...
           ~isempty(fieldnames(kalman_pll_config.(config.training_scint_model)))
            if config.is_use_cached_settings
                fprintf('Using cached values for %s.\n', config.training_scint_model);
                is_cache_used = true; % Indicate cache was used
                return; % Use cached values
            else
                fprintf('Recomputing values for %s and updating cache.\n', config.training_scint_model);
            end
        else
            fprintf(['Cache found but missing settings inside %s. ' ...
                'Computing and caching new values.\n'], config.training_scint_model);
        end
    else
        fprintf('No cache file found. Initializing kalman_pll_config.\n');
        % Initialize the kalman_pll_config struct
        kalman_pll_config = struct('CSM', struct(), 'TPPSM', struct(), 'none', struct());
        bp = 1;
    end
end

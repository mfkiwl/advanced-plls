function [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file)
% load_or_initialize_models
%
% Syntax:
%   [kalman_pll_config, is_cache_used] = load_or_initialize_models(config, cache_file)
%
% Description:
%   Loads a cached Kalman PLL configuration if available and valid.
%   If the cache does not exist or is invalid, a new configuration is initialized.
%
% Inputs:
%   config     - Struct containing user-defined PLL configuration settings.
%   cache_file - (Optional) Path to the cache file. If empty or not provided, 
%                a default path is used: 'cache/kalman_pll_cache.mat'.
%
% Outputs:
%   kalman_pll_config - Struct containing the computed or cached Kalman filter settings.
%   is_cache_used     - Boolean flag indicating whether cached values were used.
%
% Notes:
%   - If cache_file is not provided, a fixed default path is used.
%   - The function ensures that the cache directory exists before loading or saving.
%   - If cached data is invalid or missing, the configuration is reinitialized.
%
% Example:
%   config = struct(...); % user-defined configuration
%   [pll_cfg, cache_used] = load_or_initialize_models(config, '');
%   % Alternatively, specify a relative path:
%   [pll_cfg, cache_used] = load_or_initialize_models(config, 'cache/kalman_pll_cache.mat');
%
% Author:
%   Rodrigo de Lima Florindo
%   ORCID: https://orcid.org/0000-0003-0412-5583
%   Email: rdlfresearch@gmail.com

    % Default cache path if not provided
    if nargin < 2 || isempty(cache_file)
        cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');
        cache_file = fullfile(cache_dir, 'kalman_pll_cache.mat');
    else
        % Extract the directory from the provided cache file path
        cache_dir = fileparts(cache_file);
    end

    % Ensure cache directory exists
    if ~exist(cache_dir, 'dir')
        mkdir(cache_dir);
    end

    % Initialize the output flag
    is_cache_used = false;

    % Validate config input
    validateattributes(config, {'struct'}, {'nonempty'}, mfilename, 'config');

    % Load from cache if available
    if isfile(cache_file)
        fprintf('Cache file found. Loading cached kalman_pll_config.\n');
        load(cache_file, 'kalman_pll_config');

        % Check if the required configuration exists in the cache
        if isfield(kalman_pll_config, config.scintillation_training_data_config.scintillation_model) && ...
           ~isempty(fieldnames(kalman_pll_config.(config.scintillation_training_data_config.scintillation_model)))
       
            if config.is_use_cached_settings
                fprintf('Using cached values for %s.\n', config.scintillation_training_data_config.scintillation_model);
                is_cache_used = true;
                return;
            else
                fprintf('Recomputing values for %s and updating cache.\n', config.scintillation_training_data_config.scintillation_model);
            end
        else
            fprintf('Cache found but missing required settings. Computing and caching new values.\n');
        end
    else
        fprintf('No cache file found. Initializing kalman_pll_config.\n');
    end

    % Initialize a new kalman_pll_config structure
    kalman_pll_config = struct('CSM', struct(), 'TPPSM', struct(), 'none', struct());
end

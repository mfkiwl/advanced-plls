function [kalman_pll_config, is_cache_used] = get_cached_kalman_pll_config(general_config, cache_file, is_enable_cmd_print)
% get_cached_kalman_pll_general_config
%
% Syntax:
%   [kalman_pll_config, is_cache_used] = get_cached_kalman_pll_config(general_config, cache_file)
%
% Description:
%   Loads a cached Kalman PLL configuration if it exists and is valid. If the cache
%   is missing or the data is invalid, a new configuration is initialized.
%
% Inputs:
%   general_config     - A struct containing user-defined PLL configuration settings.
%   cache_file - (Optional) A string specifying the path to the cache file. If empty or
%                not provided, the default path 'cache/kalman_pll_cache.mat' is used.
%
% Outputs:
%   kalman_pll_config - A struct with the computed or loaded Kalman filter settings.
%   is_cache_used     - A boolean flag that is true if the cached configuration was used.
%   is_enable_cmd_print - Boolean flag for enabling the command prints.
%                         It is recommended to disable this option for
%                         monte carlo runs.
% Notes:
%   - If cache_file is omitted or empty, the function uses the default cache file path.
%   - The function ensures that the cache directory exists before attempting to load or
%     save the cache file.
%   - In the event that the cached data is missing or invalid, the configuration is
%     reinitialized.
%
% Example:
%   % Using the default cache file:
%   general_config = struct(...);  % User-defined configuration settings
%   [pll_cfg, cache_used] = get_cached_kalman_pll_config(general_config, '');
%
%   % Specifying a custom cache file path:
%   [pll_cfg, cache_used] = get_cached_kalman_pll_config(general_config, 'cache/kalman_pll_cache.mat');
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

    % Validate general_config input
    validateattributes(general_config, {'struct'}, {'nonempty'}, mfilename, 'general_config');

    % Load from cache if available
    if isfile(cache_file)
        if is_enable_cmd_print
            fprintf('Cache file found. Loading cached kalman_pll_config.\n');
        end
        load(cache_file, 'kalman_pll_config');

        % Check if the required configuration exists in the cache
        if isfield(kalman_pll_config, general_config.scintillation_training_data_config.scintillation_model) && ...
           ~isempty(fieldnames(kalman_pll_config.(general_config.scintillation_training_data_config.scintillation_model)))
       
            if general_config.is_use_cached_settings
                if is_enable_cmd_print
                    fprintf('Using cached values for %s.\n', general_config.scintillation_training_data_config.scintillation_model);
                end
                is_cache_used = true;
                return;
            else
                if is_enable_cmd_print
                    fprintf('Recomputing values for %s and updating cache.\n', general_config.scintillation_training_data_config.scintillation_model);
                end
            end
        else
            if is_enable_cmd_print
                fprintf('Cache found but missing required settings. Computing and caching new values.\n');
            end
        end
    else
        if is_enable_cmd_print
            fprintf('No cache file found. Initializing kalman_pll_config.\n');
        end
        
        % Initialize a new kalman_pll_config structure
        kalman_pll_config = struct('CSM', struct(), 'TPPSM', struct(), 'none', struct());
    end
end

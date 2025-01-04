% Script to analyze the output of get_csm_data and visualize results
addpath(genpath(fullfile(pwd,'..','libs')));
% Parameters for CSM data generation
S4 = 0.8;                   % Scintillation index (0 <= S4 <= 1)
tau0 = 0.7;                 % Decorrelation time (seconds)
simulation_time = 300;      % Total simulation time (seconds)
sampling_interval = 0.01;   % Sampling interval (seconds)
num_iterations = 1000;       % Number of iterations for statistical analysis

% Preallocate arrays for storing results
num_samples_per_iteration = simulation_time / sampling_interval;
psi_csm_samples = zeros(num_samples_per_iteration * num_iterations, 1);
S4_differences = zeros(1, num_iterations);

% Iterate to calculate S4 differences and collect samples
for i = 1:num_iterations
    % Generate psi_csm using the provided parameters
    psi_csm = get_csm_data(S4, tau0, simulation_time, sampling_interval);
    
    % Store the generated samples in the flattened array
    start_idx = (i - 1) * num_samples_per_iteration + 1;
    end_idx = i * num_samples_per_iteration;
    psi_csm_samples(start_idx:end_idx) = psi_csm(:);
    
    % Calculate the intensity
    intensity = abs(psi_csm).^2;
    
    % Calculate the mean of intensity and its square
    mean_intensity = mean(intensity);
    mean_intensity_squared = mean(intensity.^2);
    
    % Calculate the S4 value from the samples
    S4_calculated = sqrt((mean_intensity_squared - mean_intensity^2) / mean_intensity^2);
    
    % Store the absolute difference between calculated and input S4
    S4_differences(i) = S4_calculated - S4;
end

% Call the dedicated plot function to visualize results
plot_csm_output_analysis(S4, psi_csm_samples, S4_differences);

% Display basic statistics
mean_difference = mean(S4_differences);
std_difference = std(S4_differences);
fprintf('Mean S4 Difference: %.4f\n', mean_difference);
fprintf('Standard Deviation of S4 Differences: %.4f\n', std_difference);
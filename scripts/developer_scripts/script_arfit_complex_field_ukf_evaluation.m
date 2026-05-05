% script_arfit_complex_field_ukf_evaluation.m
%
% Developer script for evaluating the arfit-complex-field augmentation
% model with the unscented Kalman filter. The script fits a VAR model to
% the real and imaginary components of the propagated scintillation field,
% runs the UKF, and visualizes the estimated state parameters together with
% true/estimated complex-field components.
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;

script_dir = fileparts(mfilename('fullpath'));
repo_root = fullfile(script_dir, '..', '..');
addpath(genpath(fullfile(repo_root, 'libs')));

seed = 4;
rng(seed);

%% Simulation setup
doppler_profile = [0, 0, 0];
expected_doppler_profile = doppler_profile;
sampling_interval = 0.01;
L1_C_over_N0_dBHz = 42;
simulation_time = 300;
settling_time = sampling_interval;

cpssm_scenario = 'strong';
is_refractive_effects_removed_received_signal = false;
is_refractive_effects_removed_training_data = false;

[rx_sig_cpssm, los_phase, psi_cpssm] = get_received_signal( ...
    L1_C_over_N0_dBHz, ...
    'TPPSM', ...
    doppler_profile, ...
    seed, ...
    'tppsm_scenario', cpssm_scenario, ...
    'simulation_time', simulation_time, ...
    'sampling_interval', sampling_interval, ...
    'settling_time', settling_time, ...
    'is_refractive_effects_removed', is_refractive_effects_removed_received_signal);

%% UKF complex-field augmentation configuration
cache_dir = fullfile(script_dir, 'cache');
training_simulation_time = 300;
model_order = 3;
process_noise_variance_los = 1e3;

training_data_config_cpssm = struct( ...
    'scintillation_model', 'TPPSM', ...
    'scenario', cpssm_scenario, ...
    'simulation_time', training_simulation_time, ...
    'sampling_interval', sampling_interval, ...
    'is_refractive_effects_removed', is_refractive_effects_removed_training_data, ...
    'is_unwrapping_used', false);

general_config_cpssm = struct( ...
    'kf_type', 'unscented', ...
    'discrete_wiener_model_config', { {1, 3, sampling_interval, [0, 0, process_noise_variance_los], 1} }, ...
    'scintillation_training_data_config', training_data_config_cpssm, ...
    'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
    'initial_states_distributions_boundaries', { {[-pi, pi], [-0.1, 0.1], [-0.01, 0.01]} }, ...
    'expected_doppler_profile', expected_doppler_profile, ...
    'augmentation_model_initializer', struct('id', 'arfit-complex-field', 'model_params', struct('model_order', model_order)), ...
    'is_use_cached_settings', false, ...
    'is_generate_random_initial_estimates', false, ...
    'is_enable_cmd_print', false);

is_enable_cmd_print = true;
[kf_cfg, init_estimates_cpssm] = get_kalman_pll_config(general_config_cpssm, cache_dir, is_enable_cmd_print);

adaptive_cfg = struct( ...
    'measurement_cov_adapt_algorithm', 'none', ...
    'states_cov_adapt_algorithm', 'none', ...
    'sampling_interval', sampling_interval, ...
    'hard_limited', struct('is_used', false));

online_mdl_learning_cfg = struct('is_online', false);

[ukf_complex_field_cpssm, error_covariance_cpssm] = get_kalman_pll_estimates( ...
    rx_sig_cpssm, ...
    kf_cfg, ...
    init_estimates_cpssm, ...
    'unscented', ...
    'TPPSM', ...
    adaptive_cfg, ...
    online_mdl_learning_cfg);

%% Derived estimates
time_vector = sampling_interval:sampling_interval:simulation_time;
n_los_states = numel(expected_doppler_profile);
field_real_idx = n_los_states + 1;
field_imag_idx = n_los_states + 2;

field_real_est = ukf_complex_field_cpssm(:, field_real_idx);
field_imag_est = ukf_complex_field_cpssm(:, field_imag_idx);
psi_est = field_real_est + 1j * field_imag_est;
true_scint_amplitude = abs(psi_cpssm);
estimated_scint_amplitude = abs(psi_est);
true_scint_phase = get_fft_interpolated_phase(psi_cpssm);
estimated_scint_phase = get_fft_interpolated_phase(psi_est);

predicted_rx = psi_est .* exp(1j * ukf_complex_field_cpssm(:, 1));
phase_error = wrapToPi(ukf_complex_field_cpssm(:, 1) - los_phase);
field_error = psi_est - psi_cpssm;

trace_ts = zeros(size(error_covariance_cpssm, 1), 1);
for i = 1:size(error_covariance_cpssm, 1)
    trace_ts(i) = trace(squeeze(error_covariance_cpssm(i, :, :)));
end

%% Plot all estimated states
plotted_state_indices = [1:n_los_states, field_real_idx, field_imag_idx];
state_labels = make_state_labels(n_los_states);
true_los_states = get_true_los_states(time_vector, los_phase, doppler_profile, n_los_states);
true_field_states = [real(psi_cpssm), imag(psi_cpssm)];
figure('Name', 'ARFIT Complex Field UKF - State Estimates');
tiledlayout(numel(state_labels), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for plot_idx = 1:numel(plotted_state_indices)
    nexttile;
    state_idx = plotted_state_indices(plot_idx);
    hold on;
    plot(time_vector, ukf_complex_field_cpssm(:, state_idx), 'LineWidth', 1);
    if plot_idx <= size(true_los_states, 2)
        plot(time_vector, true_los_states(:, plot_idx), '--', 'LineWidth', 1);
        legend({'Estimated', 'True'}, 'Location', 'best');
    elseif state_idx == field_real_idx
        plot(time_vector, true_field_states(:, 1), '--', 'LineWidth', 1);
        legend({'Estimated', 'True'}, 'Location', 'best');
    elseif state_idx == field_imag_idx
        plot(time_vector, true_field_states(:, 2), '--', 'LineWidth', 1);
        legend({'Estimated', 'True'}, 'Location', 'best');
    end
    hold off;
    ylabel(state_labels{plot_idx}, 'Interpreter', 'none');
    grid on;
    if plot_idx == 1
        title('Estimated UKF State Parameters');
    end
    if plot_idx == numel(state_labels)
        xlabel('Time (s)');
    end
end

%% Plot true and estimated complex-field components
figure('Name', 'ARFIT Complex Field UKF - Complex Field Comparison');
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(time_vector, [real(psi_cpssm), field_real_est], 'LineWidth', 1);
ylabel('Real(\psi)');
legend({'True', 'Estimated'}, 'Location', 'best');
grid on;
title('Complex Field Components');

nexttile;
plot(time_vector, [imag(psi_cpssm), field_imag_est], 'LineWidth', 1);
ylabel('Imag(\psi)');
legend({'True', 'Estimated'}, 'Location', 'best');
grid on;

nexttile;
plot(time_vector, [real(field_error), imag(field_error)], 'LineWidth', 1);
ylabel('Field error');
legend({'Real error', 'Imag error'}, 'Location', 'best');
grid on;

nexttile;
plot(time_vector, abs(received_signal_residual(rx_sig_cpssm, predicted_rx)), 'LineWidth', 1);
ylabel('|rx error|');
xlabel('Time (s)');
grid on;

%% Plot scintillation amplitude and phase from true and estimated fields
figure('Name', 'ARFIT Complex Field UKF - Scintillation Amplitude and Phase');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(time_vector, [true_scint_amplitude, estimated_scint_amplitude], 'LineWidth', 1);
ylabel('|\psi|');
legend({'True', 'Estimated'}, 'Location', 'best');
grid on;
title('Scintillation Field Amplitude and Phase');

nexttile;
plot(time_vector, [true_scint_phase, estimated_scint_phase], 'LineWidth', 1);
ylabel('phase(\psi) (rad)');
xlabel('Time (s)');
legend({'True', 'Estimated'}, 'Location', 'best');
grid on;

%% Plot phase and covariance diagnostics
figure('Name', 'ARFIT Complex Field UKF - Diagnostics');
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;
plot(time_vector, phase_error, 'LineWidth', 1);
ylabel('LOS phase error (rad)');
grid on;

nexttile;
plot(time_vector, trace_ts, 'LineWidth', 1);
ylabel('trace(P)');
xlabel('Time (s)');
grid on;

function residual = received_signal_residual(received_signal, predicted_signal)
    residual = received_signal - predicted_signal;
end

function phase_ts = get_fft_interpolated_phase(scint_field)
    if exist('get_corrected_phase', 'file') ~= 2
        error('script_arfit_complex_field_ukf_evaluation:MissingFFTPhaseInterpolation', ...
            ['The CPSSM FFT-interpolation phase helper `get_corrected_phase` was not found on the MATLAB path. ' ...
             'Add that function to the repo/path before running this phase comparison.']);
    end
    phase_ts = get_corrected_phase(scint_field);
end

function true_los_states = get_true_los_states(time_vector, los_phase, doppler_profile, n_los_states)
    true_los_states = zeros(numel(time_vector), min(n_los_states, 3));
    if isempty(true_los_states)
        return;
    end

    true_los_states(:, 1) = los_phase;
    if size(true_los_states, 2) >= 2
        true_los_states(:, 2) = evaluate_doppler_derivative(time_vector, doppler_profile, 1);
    end
    if size(true_los_states, 2) >= 3
        true_los_states(:, 3) = evaluate_doppler_derivative(time_vector, doppler_profile, 2);
    end
end

function derivative = evaluate_doppler_derivative(time_vector, doppler_profile, derivative_order)
    derivative = zeros(numel(time_vector), 1);
    for profile_idx = derivative_order + 1:numel(doppler_profile)
        power_value = profile_idx - derivative_order - 1;
        derivative = derivative + doppler_profile(profile_idx) * ...
            (time_vector(:) .^ power_value) / factorial(power_value);
    end
end

function labels = make_state_labels(n_los_states)
    los_names = {'LOS phase', 'Doppler shift', 'Doppler drift'};
    labels = cell(1, n_los_states + 2);
    for i = 1:n_los_states
        if i <= numel(los_names)
            labels{i} = los_names{i};
        else
            labels{i} = sprintf('LOS state %d', i);
        end
    end
    labels{n_los_states + 1} = 'field real';
    labels{n_los_states + 2} = 'field imag';
end

% Script: single_freq_assessment.m
%
% Description:
%   Script for running Monte Carlo runs for performance assessement
%   of the approaches KF, AKF, KF-AR, AKF-AR, and AHL-KF-AR.
% 
% References:
%
% Author: Rodrigo de Lima Florindo
% ORCID: https://orcid.org/0000-0003-0412-5583
% Email: rdlfresearch@gmail.com

clearvars; clc;

addpath(genpath(fullfile(pwd, '..', '..', '..', 'libs')));

%inputs = get_signal_model_inputs('csm', 'strong', 1);
cache_dir = fullfile(fileparts(mfilename('fullpath')), 'cache');

[kf_ar_cfg, akf_ar_cfg, ahl_kf_ar_cfg, kf_cfg, akf_cfg, online_mdl_learning_cfg] = get_adaptive_cfgs();

% Wiener state noise variance (\sigma^2_{W,3})
sigma2_W_3_sweep = logspace(1e-14,2,5);
% Amount of Monte Carlo runs
mc_runs = 3;
% Ionospheric Scintillation Severities
severities = ["weak", "strong"];

% General parameters
% Correlation sampling interval
sampling_interval = 1e-2;
settling_time = 50;
simulation_time = 300;

% Valid time vector (after the settling périod)
valid_samples_vector = ((settling_time/sampling_interval + 1) : simulation_time/sampling_interval).';

results_matrix_template = zeros(length(sigma2_W_3_sweep),mc_runs);
struct_states = struct('phi_T', results_matrix_template, 'phi_W', results_matrix_template, 'phi_AR', results_matrix_template);
approaches_struct = struct('kf_ar', struct_states, ...
                           'akf_ar', struct_states, ...
                           'ahl_kf_ar', struct_states, ...
                           'kf', struct('phi_T',results_matrix_template), ...
                           'akf', struct('phi_T',results_matrix_template));
severities_struct = struct('weak', approaches_struct, 'strong',approaches_struct);
results = struct('csm', severities_struct, ...
    'cpssm_wo_refr', severities_struct, ...
    'cpssm_w_refr', severities_struct);

% CSM loop
for severity = severities
    for sigma2_W_3_idx = 1:length(sigma2_W_3_sweep)
        for seed = 1:mc_runs
            [rx_signal_model_inputs, gen_kf_cfg, init_estimates_csm, ~, init_estimates_none, ar_phase_idx] =...
                get_overall_cfgs(cache_dir, 'csm', severity, 'true', sigma2_W_3_sweep(sigma2_W_3_idx), sampling_interval, settling_time, simulation_time, seed);

            [rx_sig_csm, true_los_phase, ~, diffractive_phase] = get_received_signal(rx_signal_model_inputs{:});
            % For CSM, the training_scint_model is 'CSM' (AR augmented).
            % 1. KF-AR     : No adaptive update.
            % 2. AKF-AR    : nwpr adaptive update with hard_limited = false.
            % 3. AHL-KF-AR : nwpr adaptive update with hard_limited = true.
            [kf_ar_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, gen_kf_cfg, init_estimates_csm, 'standard', 'CSM', kf_ar_cfg, online_mdl_learning_cfg);
            [akf_ar_csm, ~]  = get_kalman_pll_estimates(rx_sig_csm, gen_kf_cfg, init_estimates_csm, 'standard', 'CSM', akf_ar_cfg, online_mdl_learning_cfg);
            [ahl_kf_ar_csm, ~]   = get_kalman_pll_estimates(rx_sig_csm, gen_kf_cfg, init_estimates_csm, 'standard', 'CSM', ahl_kf_ar_cfg, online_mdl_learning_cfg);
            
            % For standard KF estimates, training_scint_model is 'none'.
            % 4. KF-std    : No adaptive update.
            % 5. AKF-std   : nwpr adaptive update with hard_limited = false.
            % 6. AHL-KF-std: nwpr adaptive update with hard_limited = true.
            [kf_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, gen_kf_cfg, init_estimates_none, 'standard', 'none', kf_cfg, online_mdl_learning_cfg);
            [akf_csm, ~] = get_kalman_pll_estimates(rx_sig_csm, gen_kf_cfg, init_estimates_none, 'standard', 'none', akf_cfg, online_mdl_learning_cfg);

            % Define validation vectors (assuming valid_samples_vector is defined)
            valid_los_phase = wrapToPi(true_los_phase(valid_samples_vector,1)); % This is being iterated unecessarily.
            valid_csm_phase = diffractive_phase(valid_samples_vector,1);
            valid_total_phase = wrapToPi(valid_los_phase + valid_csm_phase);

            %%% Extracting estimates
            % KF-AR
            kf_ar_valid_hat_phi_W = kf_ar_csm(valid_samples_vector,1);
            kf_ar_valid_hat_phi_AR = kf_ar_csm(valid_samples_vector,ar_phase_idx);
            kf_ar_valid_hat_phi_T = kf_ar_valid_hat_phi_W + kf_ar_valid_hat_phi_AR;
            % AKF-AR
            akf_ar_valid_hat_phi_W = akf_ar_csm(valid_samples_vector,1);
            akf_ar_valid_hat_phi_AR = akf_ar_csm(valid_samples_vector,ar_phase_idx);
            akf_ar_valid_hat_phi_T = akf_ar_valid_hat_phi_W + akf_ar_valid_hat_phi_AR;
            % AHL-KF-AR
            ahl_kf_ar_valid_hat_phi_W = ahl_kf_ar_csm(valid_samples_vector,1);
            ahl_kf_ar_valid_hat_phi_AR = ahl_kf_ar_csm(valid_samples_vector,ar_phase_idx);
            ahl_kf_ar_valid_hat_phi_T = ahl_kf_ar_valid_hat_phi_W + ahl_kf_ar_valid_hat_phi_AR;
            % KF
            kf_valid_hat_phi_T = kf_csm(valid_samples_vector,1);
            % AKF
            akf_valid_hat_phi_T = akf_csm(valid_samples_vector,1);
            
            %%% Saving the performance assessment
            % KF-AR
            results.csm.(severity).kf_ar.phi_W(sigma2_W_3_idx, seed) = rms(wrapToPi(kf_ar_valid_hat_phi_W - valid_los_phase));
            results.csm.(severity).kf_ar.phi_AR(sigma2_W_3_idx, seed) = rms(wrapToPi(kf_ar_valid_hat_phi_AR - valid_csm_phase));
            results.csm.(severity).kf_ar.phi_T(sigma2_W_3_idx, seed) = rms(wrapToPi(kf_ar_valid_hat_phi_T - valid_total_phase));
            % AKF-AR
            results.csm.(severity).akf_ar.phi_W(sigma2_W_3_idx, seed) = rms(wrapToPi(akf_ar_valid_hat_phi_W - valid_los_phase));
            results.csm.(severity).akf_ar.phi_AR(sigma2_W_3_idx, seed) = rms(wrapToPi(akf_ar_valid_hat_phi_AR - valid_csm_phase));
            results.csm.(severity).akf_ar.phi_T(sigma2_W_3_idx, seed) = rms(wrapToPi(akf_ar_valid_hat_phi_T - valid_total_phase));
            % AHL-KF-AR
            results.csm.(severity).ahl_kf_ar.phi_W(sigma2_W_3_idx, seed) = rms(wrapToPi(ahl_kf_ar_valid_hat_phi_W - valid_los_phase));
            results.csm.(severity).ahl_kf_ar.phi_AR(sigma2_W_3_idx, seed) = rms(wrapToPi(ahl_kf_ar_valid_hat_phi_AR - valid_csm_phase));
            results.csm.(severity).ahl_kf_ar.phi_T(sigma2_W_3_idx, seed) = rms(wrapToPi(ahl_kf_ar_valid_hat_phi_T - valid_total_phase));
            % KF
            results.csm.(severity).kf.phi_T(sigma2_W_3_idx, seed) = rms(wrapToPi(kf_valid_hat_phi_T - valid_total_phase));
            % AKF
            results.csm.(severity).akf.phi_T(sigma2_W_3_idx, seed) = rms(wrapToPi(akf_valid_hat_phi_T - valid_total_phase));
        end
    end
end
% CPSSM loop (diffractive phase)

% CPSSM loop (diffractive + refractive phase)


%% Auxiliary functions
function [rx_signal_model_inputs, gen_kf_cfg, init_estimates_csm, init_estimates_cpssm, init_estimates_none, ar_phase_idx] = ...
    get_overall_cfgs(cache_dir, scint_model, severity, is_refractive_effects_removed_received_signal, sigma2_W_3, sampling_interval, settling_time, simulation_time, seed)
    % Define overall settings for the KF framework setup

    % General parameters
    doppler_profile = [0, 1000, 0.94];
    L1_C_over_N0_dBHz = 42;

    % Parameters for training the AR models for scintillation phase
    training_simulation_time = 300;
    is_refractive_effects_removed_training_data = true; % Exclude the refractive effects
    is_unwrapping_used = false; % This flag forces to use the wrapped phase for training the AR model

    % Parts for building the received signal
    csm_first_part = {L1_C_over_N0_dBHz, 'CSM', doppler_profile};
    csm_second_part = {'simulation_time', simulation_time, 'settling_time', settling_time};
    cpssm_first_part = {L1_C_over_N0_dBHz, 'TPPSM', doppler_profile};
    cpssm_second_part = {'simulation_time', simulation_time, 'settling_time', settling_time, 'is_refractive_effects_removed', is_refractive_effects_removed_received_signal};

    % CSM parameters for weak and strong scintillation
    S4_preset = [0.2, 0.9];
    tau0_preset = [1, 0.2];
    
    switch severity
        case 'weak'
            train_cfg_csm = struct('scintillation_model', 'CSM', 'S4', S4_preset(1), 'tau0', tau0_preset(1), ...
                                              'simulation_time', training_simulation_time, ...
                                              'sampling_interval', sampling_interval, ...
                                              'is_unwrapping_used', is_unwrapping_used);
            switch scint_model
                case 'csm'
                    rx_signal_model_inputs = [csm_first_part(:)',{seed},{'S4'},{S4_preset(1)},{'tau0'},{tau0_preset(1)},csm_second_part(:)'];
                    ar_model_order = 6; % See right plot of Figure 4.2 of my dissertation.
                case 'cpssm' 
                    rx_signal_model_inputs = [cpssm_first_part(:)',{seed},{'tppsm_scenario'}, {'weak'},cpssm_second_part(:)'];
                    ar_model_order = 14; % See right plot of Figure 4.7 of my dissertation.
            end
        case 'strong'
            train_cfg_csm = struct('scintillation_model', 'CSM', 'S4', S4_preset(2), 'tau0', tau0_preset(2), ...
                                              'simulation_time', training_simulation_time, ...
                                              'sampling_interval', sampling_interval, ...
                                              'is_unwrapping_used', is_unwrapping_used);
            switch scint_model
                case 'csm'
                    rx_signal_model_inputs = [csm_first_part(:)',{seed},{'S4'},{S4_preset(2)},{'tau0'},{tau0_preset(2)},csm_second_part(:)'];
                    ar_model_order = 5; % See right plot of Figure 4.2 of my dissertation.
                case 'cpssm' 
                    rx_signal_model_inputs = [cpssm_first_part(:)',{seed},{'tppsm_scenario'},{'strong'},cpssm_second_part(:)'];
                    ar_model_order = 1; % See right plot of Figure 4.7 of my dissertation.
            end
    end

    train_cfg_cpssm  = struct('scintillation_model', 'TPPSM', 'scenario', severity, ...
                              'simulation_time', training_simulation_time, ...
                              'is_refractive_effects_removed', is_refractive_effects_removed_training_data, ...
                              'sampling_interval', sampling_interval, ...
                              'is_unwrapping_used', is_unwrapping_used);
    train_cfg_none  = struct('scintillation_model', 'none', 'sampling_interval', sampling_interval);

    expected_doppler_profile = [0,1000,0.94];

    gen_cfg_csm = struct( ...
      'kf_type', 'standard', ...
      'discrete_wiener_model_config', { {1, 3, sampling_interval, [0, 0, sigma2_W_3], 1} }, ...
      'scintillation_training_data_config', train_cfg_csm, ...
      'C_over_N0_array_dBHz', L1_C_over_N0_dBHz, ...
      'initial_states_distributions_boundaries', { {[-pi, pi], [-1, 1], [-0.1, 0.1]} }, ...
      'expected_doppler_profile', expected_doppler_profile, ...
      'augmentation_model_initializer', struct('id', 'aryule', 'model_params', struct('model_order', ar_model_order)), ...
      'is_use_cached_settings', false, ...
      'is_generate_random_initial_estimates', true, ...
      'is_enable_cmd_print', false ...
    );

    gen_cfg_cpssm = gen_cfg_csm;
    gen_cfg_cpssm.scintillation_training_data_config = train_cfg_cpssm;
    
    gen_cfg_none = gen_cfg_csm;
    gen_cfg_none.scintillation_training_data_config = train_cfg_none;
    gen_cfg_none.augmentation_model_initializer.id = 'none';
    gen_cfg_none.augmentation_model_initializer.model_params = struct();
    is_enable_cmd_print = false;

    [~, init_estimates_csm] = get_kalman_pll_config(gen_cfg_csm, cache_dir, is_enable_cmd_print);
    [~, init_estimates_cpssm] = get_kalman_pll_config(gen_cfg_cpssm, cache_dir, is_enable_cmd_print);
    [gen_kf_cfg, init_estimates_none] = get_kalman_pll_config(gen_cfg_none, cache_dir, is_enable_cmd_print);

    ar_phase_idx = length(expected_doppler_profile) + 1;
end

function [kf_ar_cfg, akf_ar_cfg, ahl_kf_ar_cfg, kf_cfg, akf_cfg, online_mdl_learning_cfg] = get_adaptive_cfgs()
    % Define approaches settings
    
    sampling_interval = 1e-2;
    % Hard-Limiting constraint threshold.
    lambda = 38;  
    
    % NWPR parameters
    T_bit = 1/50;
    M_nwpr = T_bit / sampling_interval;
    N_nwpr = 20;
    
    % For AR (KFAR) estimates:
    kf_ar_cfg = struct(...
        'measurement_cov_adapt_algorithm', 'none', ...
        'states_cov_adapt_algorithm', 'none', ...
        'sampling_interval', sampling_interval, ...
        'hard_limited', struct('is_used', false));
    
    % For AR (KFAR) using the nwpr adaptive update (AKF):
    akf_ar_cfg = struct(...
        'measurement_cov_adapt_algorithm', 'nwpr', ...
        'measurement_cov_adapt_algorithm_params', struct('N_nwpr', N_nwpr, 'M_nwpr', M_nwpr), ...
        'states_cov_adapt_algorithm', 'none', ...
        'sampling_interval', sampling_interval, ...
        'hard_limited', struct('is_used', false));
    
    % For AR (KFAR) with hard limiting enabled (AHL-KF):
    ahl_kf_ar_cfg = struct(...
        'measurement_cov_adapt_algorithm', 'nwpr', ...
        'measurement_cov_adapt_algorithm_params', struct('N_nwpr', N_nwpr, 'M_nwpr', M_nwpr), ...
        'states_cov_adapt_algorithm', 'none', ...
        'sampling_interval', sampling_interval, ...
        'hard_limited', struct('is_used', true, 'L1_C_over_N0_dBHz_threshold', lambda));
    
    % For standard KF estimates (training_scint_model = 'none'):
    kf_cfg = struct(...
        'measurement_cov_adapt_algorithm', 'none', ...
        'states_cov_adapt_algorithm', 'none', ...
        'sampling_interval', sampling_interval, ...
        'hard_limited', struct('is_used', false));
    
    akf_cfg = struct(...
        'measurement_cov_adapt_algorithm', 'nwpr', ...
        'measurement_cov_adapt_algorithm_params', struct('N_nwpr', N_nwpr, 'M_nwpr', M_nwpr), ...
        'states_cov_adapt_algorithm', 'none', ...
        'sampling_interval', sampling_interval, ...
        'hard_limited', struct('is_used', false));
    
    % Online model learning setting
    online_mdl_learning_cfg = struct('is_online', false);
end
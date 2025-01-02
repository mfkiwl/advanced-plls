addpath(genpath(fullfile(pwd,'..','libs')));
C_over_N0_dBHz = 45;
S4 = 0.8;
tau0 = 0.7;
simulation_time = 300;

[psi_csm, ~] = get_received_signal(C_over_N0_dBHz,simulation_time,'MFPSM');
[psi_mfpsm, ps_realization] = get_received_signal(C_over_N0_dBHz,simulation_time,'MFPSM');
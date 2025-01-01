function [received_signal, phi_LOS] = get_received_signal(C_over_N0_dBHz,S4,tau0,simulation_time,scint_model)
% Fixed Parameters
T_I = 0.01;
phi_LOS_0 = 0;
fd = 1000;
fdr = 0.94;
rx_mean_power = 1;
B = 2e7;

phi_LOS = get_phi_LOS(simulation_time,T_I,phi_LOS_0,fd,fdr);
thermal_noise = get_thermal_noise(simulation_time,T_I,rx_mean_power, ...
    C_over_N0_dBHz,B);

switch scint_model
    case 'CSM'
        psi = get_csm_data(S4,tau0,simulation_time,T_I);
        ps_realization = NaN;
    case 'MFPSM'
        [psi, ps_realization] = get_mfpsm_data(S4,tau0,simulation_time,T_I);
    otherwise
        error(['Invalid scintillation model. Choose either `CSM` for ' ...
            'the Cornell Scintillation Model or `MFPSM` for the ' ...
            'Multi-Frequency Phase Screen Model']);
end

end
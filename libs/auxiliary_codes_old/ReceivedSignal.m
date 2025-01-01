% Project Title: EKF_VAR_Specific_Config

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code may be used to generate the specific configurations
% of the EKF-VAR for each topology.

% Its output are:
% 1 - The baseband measurements of the received signals of L1, L2 and L5 - yk;
% 2 - the true LOS phase evolution of L1, L2 and L5 - phi_D_true
% 3 - the true amplitude fades caused by ionopsheric scintillation in L1,
% L2 and L5 - rho_I_true
% 4 - the true scintillation phase of L1, L2 and L5 -  phi_I_true

function [yk,phi_D_true,rho_I_true,phi_I_true] = ReceivedSignal(TopologySelector,time,Scin_psi,AmpL1,AmpL2,AmpL5,T_I,f_receiver_Hz,fd,fdr,phiL1_0,phiL2_0,phiL5_0,seedAWGN,CN0dBHzL1,CN0dBHzL2,CN0dBHzL5,fcL1,fcL2,fcL5)
B=f_receiver_Hz/2; % Bandwith of the receiver front-end
N_int=f_receiver_Hz*T_I; % Number of samples per integration window

len = length(time);

rng(seedAWGN);

if TopologySelector == 1
    
    PbarL1 = AmpL1^2; % Power of the L1 transmitted signal
    
    %AWGN Parameters
    CN0L1 = 10^(CN0dBHzL1/10); % Signal to noise Ratio from L1 signal
    sigma2_etaL1=2*B*(PbarL1/CN0L1); 
    sigma2_eta_dL1 = sigma2_etaL1/N_int; % Variance of the AWGN noise in L1 signal
    
    phiL1_D = zeros(len,1);
    for k = 1:len
        phiL1_D(k) = phiL1_0 + 2*pi*(fd*time(k) + fdr*(time(k)^2)/2);
    end
    L1_received_signal = AmpL1*Scin_psi(:,1).*exp(1j*phiL1_D) + (randn(len,1).*sqrt(sigma2_eta_dL1/2) + 1j.*randn(len,1).*sqrt(sigma2_eta_dL1/2));
    
    yk = L1_received_signal;
    
    phi_D_true = phiL1_D;
    rho_I_true = abs(Scin_psi(:,1));
    phi_I_true = phase(Scin_psi(:,1));
    
elseif TopologySelector == 2
    
    PbarL1 = AmpL1^2; % Power of the L1 transmitted signal
    PbarL5 = AmpL5^2; % Power of the L5 transmitted signal
    
    %AWGN Parameters
    CN0L1 = 10^(CN0dBHzL1/10); % Signal to noise Ratio from L1 signal
    sigma2_etaL1=2*B*(PbarL1/CN0L1); 
    sigma2_eta_dL1 = sigma2_etaL1/N_int; % Variance of the AWGN noise in L1 signal

    CN0L5 = 10^(CN0dBHzL5/10); % Signal to noise Ratio from L2 signal
    sigma2_etaL5=2*B*(PbarL5/CN0L5); % Variance of the AWGN noise in L2 signal
    sigma2_eta_dL5 = sigma2_etaL5/N_int;
    
    phiL1_D=zeros(len,1);
    for k = 1:len
        phiL1_D(k) = phiL1_0 + 2*pi*(fd*time(k) + fdr*(time(k)^2)/2);
    end
    L1_received_signal = AmpL1*Scin_psi(:,1).*exp(1j*phiL1_D) + (randn(len,1).*sqrt(sigma2_eta_dL1/2) + 1j.*randn(len,1).*sqrt(sigma2_eta_dL1/2));

    phiL5_D=zeros(len,1);
    for k = 1:len
        phiL5_D(k) = phiL5_0 + 2*pi*(fcL5/fcL1)*(fd*time(k) + fdr*(time(k)^2)/2);
    end
    L5_received_signal = AmpL5*Scin_psi(:,2).*exp(1j*phiL5_D) + (randn(len,1).*sqrt(sigma2_eta_dL5/2) + 1j.*randn(len,1).*sqrt(sigma2_eta_dL5/2));
    
    yk = [L1_received_signal,L5_received_signal];
    
    phi_D_true = [phiL1_D,phiL5_D];
    rho_I_true = [abs(Scin_psi(:,1)),abs(Scin_psi(:,2))];
    phi_I_true = [phase(Scin_psi(:,1)),phase(Scin_psi(:,2))];
    
elseif TopologySelector == 3
    
    PbarL1 = AmpL1^2; % Power of the L1 transmitted signal
    PbarL2 = AmpL2^2; % Power of the L2 transmitted signal
    PbarL5 = AmpL5^2; % Power of the L5 transmitted signal
    
    %AWGN Parameters
    CN0L1 = 10^(CN0dBHzL1/10); % Signal to noise Ratio from L1 signal
    sigma2_etaL1=2*B*(PbarL1/CN0L1); 
    sigma2_eta_dL1 = sigma2_etaL1/N_int; % Variance of the AWGN noise in L1 signal

    CN0L2 = 10^(CN0dBHzL2/10); % Signal to noise Ratio from L2 signal
    sigma2_etaL2=2*B*(PbarL2/CN0L2); % Variance of the AWGN noise in L2 signal
    sigma2_eta_dL2 = sigma2_etaL2/N_int;

    CN0L5 = 10^(CN0dBHzL5/10); % Signal to noise Ratio from L2 signal
    sigma2_etaL5=2*B*(PbarL5/CN0L5); % Variance of the AWGN noise in L2 signal
    sigma2_eta_dL5 = sigma2_etaL5/N_int;
    
    phiL1_D=zeros(len,1);
    for k = 1:len
        phiL1_D(k) = phiL1_0 + 2*pi*(fd*time(k) + fdr*(time(k)^2)/2);
    end
    L1_received_signal = AmpL1*Scin_psi(:,1).*exp(1j*phiL1_D) + (randn(len,1).*sqrt(sigma2_eta_dL1/2) + 1j.*randn(len,1).*sqrt(sigma2_eta_dL1/2));
    
    phiL2_D=zeros(len,1);
    for k = 1:len
        phiL2_D(k) = phiL2_0 + 2*pi*(fcL2/fcL1)*(fd*time(k) + fdr*(time(k)^2)/2);
    end
    L2_received_signal = AmpL2*Scin_psi(:,2).*exp(1j*phiL2_D) + (randn(len,1).*sqrt(sigma2_eta_dL2/2) + 1j.*randn(len,1).*sqrt(sigma2_eta_dL2/2));
    
    phiL5_D=zeros(len,1);
    for k = 1:len
        phiL5_D(k) = phiL5_0 + 2*pi*(fcL5/fcL1)*(fd*time(k) + fdr*(time(k)^2)/2);
    end
    L5_received_signal = AmpL5*Scin_psi(:,3).*exp(1j*phiL5_D) + (randn(len,1).*sqrt(sigma2_eta_dL5/2) + 1j.*randn(len,1).*sqrt(sigma2_eta_dL5/2));
    
    yk = [L1_received_signal,L2_received_signal,L5_received_signal];
    
    phi_D_true = [phiL1_D,phiL2_D,phiL5_D];
    rho_I_true = [abs(Scin_psi(:,1)),abs(Scin_psi(:,2)),abs(Scin_psi(:,3))];
    phi_I_true = [phase(Scin_psi(:,1)),phase(Scin_psi(:,2)),phase(Scin_psi(:,3))];
end
end
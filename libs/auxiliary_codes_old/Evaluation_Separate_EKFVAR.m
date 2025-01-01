% Project Title: Evaulation_EKFVAR

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 03/07/2024

% Description: This code may be used to evaluate all states of interest
% estimated by the EKF-VAR for each of its available topologies.

% Its output are:
% 1 - The time series of the RMSE of the LOS phases - RMSED
% 1 - The time series of the RMSE of the scintillation phases - RMSEI
% 1 - The time series of the RMSE of the total carrier phases - RMSET
% 1 - The time series of the RMSE of the scintillation amplitudes - RMSERho
% 1 - The time series of the residuals of the LOS phases - errorD
% 1 - The time series of the residuals of the scintillation phases - errorI
% 1 - The time series of the residuals of the total carrier phases - errorT
% 1 - The time series of the residuals of the scintillation amplitudes - errorRho

function [RMSED,RMSEI,RMSET,RMSERho,errorD,errorI,errorT,errorRho] = Evaluation_Separate_EKFVAR(xkk1,phi_D_true,phi_I_true,rho_I_true,TopologySelector,L,M,P_rho,P_phi)
    if TopologySelector == 1
        phiL1_D_error = xkk1(:,1) - phi_D_true(:,1);
        phiL1_T_error = xkk1(:,1) + xkk1(:,(L+M)+L*P_rho) - phi_D_true(:,1) - phi_I_true(:,1);
        phiL1_I_error = phiL1_T_error - phiL1_D_error;
        rhoL1_error = xkk1(:,L+M) - rho_I_true(:,1);
        
        RMSE_L1_D = sqrt(mean(phiL1_D_error.^2));
        RMSE_L1_I = sqrt(mean(phiL1_I_error.^2));
        RMSE_L1_T = sqrt(mean(phiL1_T_error.^2));
        RMSE_L1_rho = sqrt(mean(rhoL1_error.^2));
        
        %Outputs
        errorD = phiL1_D_error;
        errorI = phiL1_I_error;
        errorT = phiL1_T_error;
        errorRho = rhoL1_error;
        
        RMSED = RMSE_L1_D;
        RMSEI = RMSE_L1_I;
        RMSET = RMSE_L1_T;
        RMSERho = RMSE_L1_rho;
    elseif TopologySelector == 2
        phiL1_D_error = xkk1(:,1) - phi_D_true(:,1);
        phiL1_T_error = xkk1(:,1) + xkk1(:,(L+M) + L*P_rho) - phi_D_true(:,1) - phi_I_true(:,1);
        phiL1_I_error = phiL1_T_error - phiL1_D_error;
        rhoL1_error = xkk1(:,L+M) - rho_I_true(:,1);
        
        RMSE_L1_D = sqrt(mean(phiL1_D_error.^2));
        RMSE_L1_I = sqrt(mean(phiL1_I_error.^2));
        RMSE_L1_T = sqrt(mean(phiL1_T_error.^2));
        RMSE_L1_rho = sqrt(mean(rhoL1_error.^2));

        phiL5_D_error = xkk1(:,2) - phi_D_true(:,2);
        phiL5_T_error = xkk1(:,2) + xkk1(:,(L+M) + L*P_rho + 1) - phi_D_true(:,2) - phi_I_true(:,2);
        phiL5_I_error = phiL5_T_error - phiL5_D_error;
        rhoL5_error = xkk1(:,L+M+1) - rho_I_true(:,2);

        RMSE_L5_D = sqrt(mean(phiL5_D_error.^2));
        RMSE_L5_I = sqrt(mean(phiL5_I_error.^2));
        RMSE_L5_T = sqrt(mean(phiL5_T_error.^2));
        RMSE_L5_rho = sqrt(mean(rhoL5_error.^2));
        
        %Outputs
        errorD = [phiL1_D_error,phiL5_D_error];
        errorI = [phiL1_I_error,phiL5_I_error];
        errorT = [phiL1_T_error,phiL5_T_error];
        errorRho = [rhoL1_error,rhoL5_error];
        
        RMSED = [RMSE_L1_D,RMSE_L5_D];
        RMSEI = [RMSE_L1_I,RMSE_L5_I];
        RMSET = [RMSE_L1_T,RMSE_L5_T];
        RMSERho = [RMSE_L1_rho,RMSE_L5_rho];
    elseif TopologySelector == 3
        phiL1_D_error = xkk1(:,1) - phi_D_true(:,1);
        phiL1_T_error = xkk1(:,1) + xkk1(:,(L+M)+L*P_rho) - phi_D_true(:,1) - phi_I_true(:,1);
        phiL1_I_error = phiL1_T_error - phiL1_D_error;
        rhoL1_error = xkk1(:,L+M) - rho_I_true(:,1);
        
        RMSE_L1_D = sqrt(mean(phiL1_D_error.^2));
        RMSE_L1_I = sqrt(mean(phiL1_I_error.^2));
        RMSE_L1_T = sqrt(mean(phiL1_T_error.^2));
        RMSE_L1_rho = sqrt(mean(rhoL1_error.^2));

        phiL2_D_error = xkk1(:,2) - phi_D_true(:,2);
        phiL2_T_error = xkk1(:,2) + xkk1(:,(L+M) + L*P_rho + 1) - phi_D_true(:,2) - phi_I_true(:,2);
        phiL2_I_error = phiL2_T_error - phiL2_D_error;
        rhoL2_error = xkk1(:,L+M+1) - rho_I_true(:,2);

        RMSE_L2_D = sqrt(mean(phiL2_D_error.^2));
        RMSE_L2_I = sqrt(mean(phiL2_I_error.^2));
        RMSE_L2_T = sqrt(mean(phiL2_T_error.^2));
        RMSE_L2_rho = sqrt(mean(rhoL2_error.^2));

        phiL5_D_error = xkk1(:,3) - phi_D_true(:,3);
        phiL5_T_error = xkk1(:,3) + xkk1(:,(L+M) + L*P_rho + 2) - phi_D_true(:,3) - phi_I_true(:,3);
        phiL5_I_error = phiL5_T_error - phiL5_D_error;
        rhoL5_error = xkk1(:,L+M+2) - rho_I_true(:,3);
        
        RMSE_L5_D = sqrt(mean(phiL5_D_error.^2));
        RMSE_L5_I = sqrt(mean(phiL5_I_error.^2));
        RMSE_L5_T = sqrt(mean(phiL5_T_error.^2));
        RMSE_L5_rho = sqrt(mean(rhoL5_error.^2));
        
        %Outputs
        errorD = [phiL1_D_error,phiL2_D_error,phiL5_D_error];
        errorI = [phiL1_I_error,phiL2_I_error,phiL5_I_error];
        errorT = [phiL1_T_error,phiL2_T_error,phiL5_T_error];
        errorRho = [rhoL1_error,rhoL2_error,rhoL5_error];
        
        RMSED = [RMSE_L1_D,RMSE_L2_D,RMSE_L5_D];
        RMSEI = [RMSE_L1_I,RMSE_L2_I,RMSE_L5_I];
        RMSET = [RMSE_L1_T,RMSE_L2_T,RMSE_L5_T];
        RMSERho = [RMSE_L1_rho,RMSE_L2_rho,RMSE_L5_rho];
    end
end
% Project Title: Evaulation_KFVAR

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code may be used to evaluate all states of interest
% estimated by the KF-VAR for each of its available topologies.

% Its output are:
% 1 - The time series of the RMSE of the LOS phases - RMSED
% 1 - The time series of the RMSE of the scintillation phases - RMSEI
% 1 - The time series of the RMSE of the total carrier phases - RMSET
% 1 - The time series of the residuals of the LOS phases - errorD
% 1 - The time series of the residuals of the scintillation phases - errorI
% 1 - The time series of the residuals of the total carrier phases - errorT

function [RMSED,RMSEI,RMSET,errorD,errorI,errorT] = Evaluation_KFVAR(xkk1,phi_D_true,phi_I_true,TopologySelector,L,M)
    if TopologySelector == 1
        phiL1_D_error = xkk1(:,1) - phi_D_true(:,1);
        phiL1_T_error = xkk1(:,1) + xkk1(:,L+M) - phi_D_true(:,1) - phi_I_true(:,1);
        phiL1_I_error = phiL1_T_error - phiL1_D_error;

        RMSE_L1_D = sqrt(mean(phiL1_D_error.^2));
        RMSE_L1_I = sqrt(mean(phiL1_I_error.^2));
        RMSE_L1_T = sqrt(mean(phiL1_T_error.^2));
        
        %Outputs
        errorD = phiL1_D_error;
        errorI = phiL1_I_error;
        errorT = phiL1_T_error;
        
        RMSED = RMSE_L1_D;
        RMSEI = RMSE_L1_I;
        RMSET = RMSE_L1_T;
    elseif TopologySelector == 2
        phiL1_D_error = xkk1(:,1) - phi_D_true(:,1);
        phiL1_T_error = xkk1(:,1) + xkk1(:,L+M) - phi_D_true(:,1) - phi_I_true(:,1);
        phiL1_I_error = phiL1_T_error - phiL1_D_error;

        RMSE_L1_D = sqrt(mean(phiL1_D_error.^2));
        RMSE_L1_I = sqrt(mean(phiL1_I_error.^2));
        RMSE_L1_T = sqrt(mean(phiL1_T_error.^2));

        phiL5_D_error = xkk1(:,2) - phi_D_true(:,2);
        phiL5_T_error = xkk1(:,2) + xkk1(:,L+M+1) - phi_D_true(:,2) - phi_I_true(:,2);
        phiL5_I_error = phiL5_T_error - phiL5_D_error;

        RMSE_L5_D = sqrt(mean(phiL5_D_error.^2));
        RMSE_L5_I = sqrt(mean(phiL5_I_error.^2));
        RMSE_L5_T = sqrt(mean(phiL5_T_error.^2));
        
        %Outputs
        errorD = [phiL1_D_error,phiL5_D_error];
        errorI = [phiL1_I_error,phiL5_I_error];
        errorT = [phiL1_T_error,phiL5_T_error];
        
        RMSED = [RMSE_L1_D,RMSE_L5_D];
        RMSEI = [RMSE_L1_I,RMSE_L5_I];
        RMSET = [RMSE_L1_T,RMSE_L5_T];
    elseif TopologySelector == 3
        phiL1_D_error = xkk1(:,1) - phi_D_true(:,1);
        phiL1_T_error = xkk1(:,1) + xkk1(:,L+M) - phi_D_true(:,1) - phi_I_true(:,1);
        phiL1_I_error = phiL1_T_error - phiL1_D_error;

        RMSE_L1_D = sqrt(mean(phiL1_D_error.^2));
        RMSE_L1_I = sqrt(mean(phiL1_I_error.^2));
        RMSE_L1_T = sqrt(mean(phiL1_T_error.^2));

        phiL2_D_error = xkk1(:,2) - phi_D_true(:,2);
        phiL2_T_error = xkk1(:,2) + xkk1(:,L+M+1) - phi_D_true(:,2) - phi_I_true(:,2);
        phiL2_I_error = phiL2_T_error - phiL2_D_error;

        RMSE_L2_D = sqrt(mean(phiL2_D_error.^2));
        RMSE_L2_I = sqrt(mean(phiL2_I_error.^2));
        RMSE_L2_T = sqrt(mean(phiL2_T_error.^2));

        phiL5_D_error = xkk1(:,3) - phi_D_true(:,3);
        phiL5_T_error = xkk1(:,3) + xkk1(:,L+M+2) - phi_D_true(:,3) - phi_I_true(:,3);
        phiL5_I_error = phiL5_T_error - phiL5_D_error;

        RMSE_L5_D = sqrt(mean(phiL5_D_error.^2));
        RMSE_L5_I = sqrt(mean(phiL5_I_error.^2));
        RMSE_L5_T = sqrt(mean(phiL5_T_error.^2));
        
        %Outputs
        errorD = [phiL1_D_error,phiL2_D_error,phiL5_D_error];
        errorI = [phiL1_I_error,phiL2_I_error,phiL5_I_error];
        errorT = [phiL1_T_error,phiL2_T_error,phiL5_T_error];
        
        RMSED = [RMSE_L1_D,RMSE_L2_D,RMSE_L5_D];
        RMSEI = [RMSE_L1_I,RMSE_L2_I,RMSE_L5_I];
        RMSET = [RMSE_L1_T,RMSE_L2_T,RMSE_L5_T];
    end
end
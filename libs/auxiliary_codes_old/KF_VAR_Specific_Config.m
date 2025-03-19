% Project Title: KF_VAR_Specific_Config

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code may be used to generate the specific configurations
% of the KF-VAR for each topology.

% Its output are:
% 1 - A matrix containing the complex values related to the ionospheric
% scintillation realizations of MFPSM - Scin_psi;
% 2 - the bias vector of the VAR model - c;
% 3 - the regressor parameters of the VAR model - A;
% 4 - the covariance matrix of the VAR model - Q_arfit;
% 5 - the number of frequencies being tracked - L;
% 6 - the parameter of relations between the frequencies that is used to
% compute the compute the state transition sub-matrix of the LOS dynamics
% state space - deltaArray
% 7 - the number of states from the VAR model - NumSeries;
% 8 - the VAR model order - P;
% 9 - the measurement noise covariance matrix - Rk;
% 10 - Auxiliary arrays used to select the phase
% scintillation estimates states (Aux_thetaL1,Aux_thetaL2,Aux_thetaL5)
% 11 - The measurement state transition matrix - Hk

if TopologySelector == 1
    Scin_psi(:,1) = Y_obs_full(:,1,scintSeed).*exp(1j * Y_obs_full(:,4,scintSeed));
    [c,A,Q_arfit,~,~,~] = arfit(Y_obs_full(:,4,scintSeed),1,30);
    L = 1;
    deltaArray = 1;
    
    NumSeries = size(A,1); % Number of states from the VAR model 
    P = size(A,2)/NumSeries; % VAR model order
    
    Rk = sigma2Measure_valueL1;
    
    Aux_thetaL1=[1,zeros(1,M-1),1,zeros(1,(NumSeries*P-1))];
    Hk = Aux_thetaL1;
elseif TopologySelector == 2
    Scin_psi(:,1) = Y_obs_full(:,1,scintSeed).*exp(1j * Y_obs_full(:,4,scintSeed));
    Scin_psi(:,2) = Y_obs_full(:,3,scintSeed).*exp(1j * Y_obs_full(:,6,scintSeed));
    [c,A,Q_arfit,~,~,~] = arfit(Y_obs_full(:,[4,6],scintSeed),1,30);
    L = 2;
    deltaArray = [1,fcL5/fcL1];
    
    NumSeries = size(A,1); % Number of states from the VAR model
    P = size(A,2)/NumSeries; % VAR model order
    
    Rk = diag([sigma2Measure_valueL1,sigma2Measure_valueL5]);
    
    Aux_thetaL1=[1,0,zeros(1,M-1),1,zeros(1,(NumSeries*P-1))];
    Aux_thetaL5=[0,1,zeros(1,M-1),0,1,zeros(1,(NumSeries*P-2))];
    Hk = [Aux_thetaL1;Aux_thetaL5];
elseif TopologySelector == 3
    Scin_psi(:,1) = Y_obs_full(:,1,scintSeed).*exp(1j * Y_obs_full(:,4,scintSeed));
    Scin_psi(:,2) = Y_obs_full(:,2,scintSeed).*exp(1j * Y_obs_full(:,5,scintSeed));
    Scin_psi(:,3) = Y_obs_full(:,3,scintSeed).*exp(1j * Y_obs_full(:,6,scintSeed));
    [c,A,Q_arfit,~,~,~] = arfit(Y_obs_full(:,4:6,scintSeed),1,30);
    L = 3;
    deltaArray = [1,fcL2/fcL1,fcL5/fcL1];
    
    NumSeries = size(A,1); % Number of states from the VAR model 
    P = size(A,2)/NumSeries; % VAR model order
    
    Rk = diag([sigma2Measure_valueL1,sigma2Measure_valueL2,sigma2Measure_valueL5]);
    
    Aux_thetaL1=[1,0,0,zeros(1,M-1),1,zeros(1,(NumSeries*P-1))];
    Aux_thetaL2=[0,1,0,zeros(1,M-1),0,1,zeros(1,(NumSeries*P-2))];
    Aux_thetaL5=[0,0,1,zeros(1,M-1),0,0,1,zeros(1,(NumSeries*P-3))];
    Hk = [Aux_thetaL1;Aux_thetaL2;Aux_thetaL5];
end
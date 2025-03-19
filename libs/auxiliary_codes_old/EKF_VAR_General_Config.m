%% Comments
% Project Title: EKF_VAR_General_Config

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code may be used to generate the initial configurations
% of the EKF-VAR for each topology.

% Its output are:
% 1 - The full initial error covariance matrix - Px_0_0; 
% 2 - The initial States - x_hat_0_0;
% 3 - the full process noise covariance matrix - Qk;
% 4 - the full process state transition matrix - Fk;
% 5 - the Constant of the VAR model - C.

function [Px_0_0,x_hat_0_0,Qk,Fk,C] = EKF_VAR_General_Config(L,M,NumSeries,P,T_I,sigma2,fd,fdr,deltaArray,Q_arfit,A,c,TopologySelector)
sigma_QD = zeros(1,L+M-1);
sigma_QD(end) = sigma2;%2*pi*1.35*10^-8;
[FD,QD]=LOS(L,M,T_I,sigma_QD,deltaArray);

%Initialization Matrices
DopplerInit_States = zeros(1,(L+M-1));
DopplerInit_States(L+1) = fd;
DopplerInit_States(L+2) = fdr;
DopplerInit_Covariance = zeros(1,(L+M-1));

if TopologySelector == 1
    DopplerInit_Covariance(1) = pi^2/3;
    DopplerInit_Covariance(2) = 5^2/12;
    DopplerInit_Covariance(3) = 0.01^2/12;
elseif TopologySelector == 2
    DopplerInit_Covariance(1) = pi^2/3;
    DopplerInit_Covariance(2) = pi^2/3;
    DopplerInit_Covariance(3) = 5^2/12;
    DopplerInit_Covariance(4) = 0.01^2/12;
elseif TopologySelector == 3
    DopplerInit_Covariance(1) = pi^2/3;
    DopplerInit_Covariance(2) = pi^2/3;
    DopplerInit_Covariance(3) = pi^2/3;
    DopplerInit_Covariance(4) = 5^2/12;
    DopplerInit_Covariance(5) = 0.01^2/12;
end

ScintInit_States = [ones(L,P);zeros(L,P)];
ScintInit_States = ScintInit_States(:).';
ScintInit_Covariance = [ones(L,P)*((0.1^2)/12);ones(L,P)*(pi^2/3)];
ScintInit_Covariance = ScintInit_Covariance(:).';

Px_0_0 = diag([DopplerInit_Covariance,ScintInit_Covariance]); % Initial Error Covariance Matrix
x_hat_0_0 = ([DopplerInit_States,ScintInit_States])'; % Initial state vector

QVAR = blkdiag(Q_arfit,diag(zeros(1,NumSeries*(P-1))));
Qk = blkdiag(QD,QVAR);

FVAR = [A;[eye(NumSeries*(P-1)),zeros(NumSeries*(P-1),NumSeries)]];

Fk = blkdiag(FD,FVAR);

C = [zeros(1,L+M-1),c',zeros(1, NumSeries*(P-1))];
end
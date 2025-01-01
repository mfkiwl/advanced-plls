% Project Title: KF_VAR_PLL

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code run a KF algorithm for the proposed SSM.

% Its output are:
% 1 - The time series of the estimated states using knowledge until instant k-1 - xkk1_TS
% 2 - The time series of the estimated error covariance matrix using knowledge until instant k-1 - Pkk1_TS;

function [xkk1_TS, Pkk1_TS] = KF_VAR_PLL(Px_0_0,x_hat_0_0,yk,Fk,Qk,C,Rk,Hk,AdaptSwitch,T_I,CN0L1,CN0L2,CN0L5,TopologySelector)
Samples = 30001;
xkk1_TS = zeros(Samples,length(x_hat_0_0));
Pkk1_TS = zeros(size(Px_0_0,1),size(Px_0_0,2),Samples);

for k = 1:Samples
    %% Prediction
    
    if k == 1 % Initialization
        xkk1 = x_hat_0_0;
        Pkk1 = Px_0_0;
    else
        xk1k1 = xkk;
        Pk1k1 = Pkk;
        
        xkk1 = Fk*xk1k1 + C'; % State propagation
        Pkk1 = Fk*Pk1k1*Fk' + Qk; % Error covariance propagation
    end
    xkk1_TS(k,:) = xkk1';
    Pkk1_TS(:,:,k) = Pkk1;
    %% Atualization
    ykk1_phi_total = Hk*xkk1;
    if TopologySelector == 1
        ykk1 = exp(-1j*ykk1_phi_total(1));
    elseif TopologySelector == 2
        ykk1 = [exp(-1j*ykk1_phi_total(1));exp(-1j*ykk1_phi_total(2))];
    elseif TopologySelector == 3
        ykk1 = [exp(-1j*ykk1_phi_total(1));exp(-1j*ykk1_phi_total(2));exp(-1j*ykk1_phi_total(3))];
    end
    
    if AdaptSwitch == 1
        Rk_adapt = Rk;
    elseif AdaptSwitch == 2
        if TopologySelector == 1
            CN0_L1_hat = CN0L1*abs(yk(k,1))^2;
            sigma2Measure_valueL1 = (1/(2*CN0_L1_hat*T_I))*(1+(1/(2*CN0_L1_hat*T_I)));
            Rk_adapt = sigma2Measure_valueL1;
        elseif TopologySelector == 2
            CN0_L1_hat = CN0L1*abs(yk(k,1))^2;
            CN0_L5_hat = CN0L5*abs(yk(k,2))^2;
            sigma2Measure_valueL1 = (1/(2*CN0_L1_hat*T_I))*(1+(1/(2*CN0_L1_hat*T_I)));
            sigma2Measure_valueL5 = (1/(2*CN0_L5_hat*T_I))*(1+(1/(2*CN0_L5_hat*T_I)));
            Rk_adapt = diag([sigma2Measure_valueL1,sigma2Measure_valueL5]);
        elseif TopologySelector == 3
            CN0_L1_hat = CN0L1*abs(yk(k,1))^2;
            CN0_L2_hat = CN0L2*abs(yk(k,2))^2;
            CN0_L5_hat = CN0L5*abs(yk(k,3))^2;
            kp = 1;
            sigma2Measure_valueL1 = kp*(1/(2*CN0_L1_hat*T_I))*(1+(1/(2*CN0_L1_hat*T_I)));
            sigma2Measure_valueL2 = kp*(1/(2*CN0_L2_hat*T_I))*(1+(1/(2*CN0_L2_hat*T_I)));
            sigma2Measure_valueL5 = kp*(1/(2*CN0_L5_hat*T_I))*(1+(1/(2*CN0_L5_hat*T_I)));
            Rk_adapt = diag([sigma2Measure_valueL1,sigma2Measure_valueL2,sigma2Measure_valueL5]);
        end
    end

    Kk = Pkk1*Hk'*((Hk*Pkk1*Hk' + Rk_adapt)\eye(TopologySelector)); % Kalman Gain computation

    if TopologySelector == 1
        ik_aux = zeros(1,1);
        ik_aux(1,1) = yk(k,1)*ykk1(1);
        ik = atan2(imag(ik_aux(1)),real(ik_aux(1)));
    elseif TopologySelector == 2
        ik_aux = zeros(2,1);
        ik_aux(1,1) = yk(k,1)*ykk1(1);
        ik_aux(2,1) = yk(k,2)*ykk1(2);
        ik = [atan2(imag(ik_aux(1)),real(ik_aux(1)));atan2(imag(ik_aux(2)),real(ik_aux(2)))];
    elseif TopologySelector == 3
        ik_aux = zeros(3,1);
        ik_aux(1,1) = yk(k,1)*ykk1(1);
        ik_aux(2,1) = yk(k,2)*ykk1(2);
        ik_aux(3,1) = yk(k,3)*ykk1(3);
        ik = [atan2(imag(ik_aux(1)),real(ik_aux(1)));atan2(imag(ik_aux(2)),real(ik_aux(2)));atan2(imag(ik_aux(3)),real(ik_aux(3)))];
    end

    xkk = xkk1 + Kk*ik;
    Pkk = Pkk1 - Kk*Hk*Pkk1;
end
end
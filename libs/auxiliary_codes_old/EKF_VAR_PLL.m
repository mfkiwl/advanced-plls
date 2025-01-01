% Project Title: EKF_VAR_PLL

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code run a EKF algorithm for the proposed SSM.

% Its output are:
% 1 - The time series of the estimated states using knowledge until instant k-1 - xkk1_TS
% 2 - The time series of the estimated error covariance matrix using knowledge until instant k-1 - Pkk1_TS;

function [xkk1_TS, Pkk1_TS] = EKF_VAR_PLL(Px_0_0,x_hat_0_0,yk,Rk,Fk,Qk,M,L,NumSeries,P,AmpL1,AmpL2,AmpL5,Aux,C,TopologySelector,AdaptSwitch,CN0L1,CN0L2,CN0L5,T_I)

Samples = 30001;

xkk1_TS = zeros(Samples,length(x_hat_0_0));
Pkk1_TS = zeros(size(Px_0_0,1),size(Px_0_0,2),Samples);

yk_formatted = zeros(Samples,NumSeries);

if TopologySelector == 1
    yk_formatted(:,1) = real(yk(:,1));
    yk_formatted(:,2) = imag(yk(:,1));
elseif TopologySelector == 2
    yk_formatted(:,1) = real(yk(:,1));
    yk_formatted(:,2) = imag(yk(:,1));
    yk_formatted(:,3) = real(yk(:,2));
    yk_formatted(:,4) = imag(yk(:,2));
elseif TopologySelector == 3
    yk_formatted(:,1) = real(yk(:,1));
    yk_formatted(:,2) = imag(yk(:,1));
    yk_formatted(:,3) = real(yk(:,2));
    yk_formatted(:,4) = imag(yk(:,2));
    yk_formatted(:,5) = real(yk(:,3));
    yk_formatted(:,6) = imag(yk(:,3));
end

for k = 1:Samples
    %% Prediction
    if k == 1
        xkk1 = x_hat_0_0;
        Pkk1 = Px_0_0;
    else
        xk1k1 = xkk;
        Pk1k1 = Pkk;
        xkk1 = Fk*xk1k1 + C';
        Pkk1 = Fk*Pk1k1*Fk' + Qk;
    end

    % There is a trade-off in choosing higher thresholds for the amplitude
    % fades. When gamma is higher than 0.01, the Total carrier phase RMSE
    % closer to the KF, but the LOS phase and Scint Phase are worse. When
    % it is below 0.01, it gets
    
    % Amplitude Constraints
    gamma = 0.1;
    if TopologySelector == 1
        if xkk1(L+M-1+1) <= gamma
            xkk1(L+M-1+1) = gamma;
        end
    elseif TopologySelector == 2
        if xkk1(L+M-1+1) <= gamma
            xkk1(L+M-1+1) = gamma;
        end
        if xkk1(L+M-1+2) <= gamma
            xkk1(L+M-1+2) = gamma;
        end
    elseif TopologySelector == 3
        if xkk1(L+M-1+1) <= gamma
            xkk1(L+M-1+1) = gamma;
        end
        if xkk1(L+M-1+2) <= gamma
            xkk1(L+M-1+2) = gamma;
        end
        if xkk1(L+M-1+3) <= gamma
            xkk1(L+M-1+3) = gamma;
        end
    end
    
    xkk1_TS(k,:) = xkk1';
    Pkk1_TS(:,:,k) = Pkk1;
    %% Atualization

    % Jacobian Calculation
    JHd = zeros(NumSeries,L);
    JHrho = zeros(NumSeries,L);
    JHphi = zeros(NumSeries,L);
    ykk1_aux = Aux*xkk1;
    if TopologySelector == 1
        JHd(1,1) = -AmpL1*ykk1_aux(1)*sin(ykk1_aux(2));
        JHd(2,1) = AmpL1*ykk1_aux(1)*cos(ykk1_aux(2));

        JHrho(1,1) = AmpL1*cos(ykk1_aux(2));
        JHrho(2,1) = AmpL1*sin(ykk1_aux(2));

        JHphi = JHd;
        
        I_L1_kk1 = real(ykk1_aux(1)*AmpL1*exp(1j*ykk1_aux(2)));
        Q_L1_kk1 = imag(ykk1_aux(1)*AmpL1*exp(1j*ykk1_aux(2)));

        ykk1 = [I_L1_kk1;Q_L1_kk1];
    elseif TopologySelector == 2

        JHd(1,1) = -AmpL1*ykk1_aux(1)*sin(ykk1_aux(3));
        JHd(2,1) = AmpL1*ykk1_aux(1)*cos(ykk1_aux(3));
        JHd(3,2) = -AmpL5*ykk1_aux(2)*sin(ykk1_aux(4));
        JHd(4,2) = AmpL5*ykk1_aux(2)*cos(ykk1_aux(4));


        JHrho(1,1) = AmpL1*cos(ykk1_aux(3));
        JHrho(2,1) = AmpL1*sin(ykk1_aux(3));
        JHrho(3,2) = AmpL5*cos(ykk1_aux(4));
        JHrho(4,2) = AmpL5*sin(ykk1_aux(4));

        JHphi = JHd;

        I_L1_kk1 = real(ykk1_aux(1)*exp(1j*ykk1_aux(3)));
        Q_L1_kk1 = imag(ykk1_aux(1)*exp(1j*ykk1_aux(3)));
        I_L5_kk1 = real(ykk1_aux(2)*exp(1j*ykk1_aux(4)));
        Q_L5_kk1 = imag(ykk1_aux(2)*exp(1j*ykk1_aux(4)));

        ykk1 = [I_L1_kk1;Q_L1_kk1;I_L5_kk1;Q_L5_kk1];
    elseif TopologySelector == 3
        
        JHd(1,1) = -AmpL1*ykk1_aux(1)*sin(ykk1_aux(4));
        JHd(2,1) = AmpL1*ykk1_aux(1)*cos(ykk1_aux(4));
        JHd(3,2) = -AmpL2*ykk1_aux(2)*sin(ykk1_aux(5));
        JHd(4,2) = AmpL2*ykk1_aux(2)*cos(ykk1_aux(5));
        JHd(5,3) = -AmpL5*ykk1_aux(3)*sin(ykk1_aux(6));
        JHd(6,3) = AmpL5*ykk1_aux(3)*cos(ykk1_aux(6));
        
        JHrho(1,1) = AmpL1*cos(ykk1_aux(4));
        JHrho(2,1) = AmpL1*sin(ykk1_aux(4));
        JHrho(3,2) = AmpL2*cos(ykk1_aux(5));
        JHrho(4,2) = AmpL2*sin(ykk1_aux(5));
        JHrho(5,3) = AmpL5*cos(ykk1_aux(6));
        JHrho(6,3) = AmpL5*sin(ykk1_aux(6));
        
        JHphi = JHd;
        
        I_L1_kk1 = real(ykk1_aux(1)*exp(1j*ykk1_aux(4)));
        Q_L1_kk1 = imag(ykk1_aux(1)*exp(1j*ykk1_aux(4)));
        I_L2_kk1 = real(ykk1_aux(2)*exp(1j*ykk1_aux(5)));
        Q_L2_kk1 = imag(ykk1_aux(2)*exp(1j*ykk1_aux(5)));
        I_L5_kk1 = real(ykk1_aux(3)*exp(1j*ykk1_aux(6)));
        Q_L5_kk1 = imag(ykk1_aux(3)*exp(1j*ykk1_aux(6)));
        
        ykk1 = [I_L1_kk1;Q_L1_kk1;I_L2_kk1;Q_L2_kk1;I_L5_kk1;Q_L5_kk1];
    end
    
    JHk = [JHd,zeros(NumSeries,M-1),JHrho,JHphi,zeros(NumSeries,NumSeries*(P-1))];
    
    if AdaptSwitch == 1
        Rk_adapt = Rk;
    elseif AdaptSwitch == 2
        if TopologySelector == 1
            CN0_L1_hat = CN0L1*abs(yk(1))^2;
            sigma2Measure_valueL1 = (1/(2*CN0_L1_hat*T_I))*(1+(1/(2*CN0_L1_hat*T_I)));
            Rk_adapt = diag([sigma2Measure_valueL1/2,sigma2Measure_valueL1/2]);
        elseif TopologySelector == 2
            CN0_L1_hat = CN0L1*abs(yk(1))^2;
            CN0_L5_hat = CN0L5*abs(yk(2))^2;
            sigma2Measure_valueL1 = (1/(2*CN0_L1_hat*T_I))*(1+(1/(2*CN0_L1_hat*T_I)));
            sigma2Measure_valueL5 = (1/(2*CN0_L5_hat*T_I))*(1+(1/(2*CN0_L5_hat*T_I)));
            Rk_adapt = diag([sigma2Measure_valueL1/2,sigma2Measure_valueL1/2,sigma2Measure_valueL5/2,sigma2Measure_valueL5/2]);
        elseif TopologySelector == 3
            CN0_L1_hat = CN0L1*abs(yk(1))^2;
            CN0_L2_hat = CN0L2*abs(yk(2))^2;
            CN0_L5_hat = CN0L5*abs(yk(3))^2;
            sigma2Measure_valueL1 = (1/(2*CN0_L1_hat*T_I))*(1+(1/(2*CN0_L1_hat*T_I)));
            sigma2Measure_valueL2 = (1/(2*CN0_L2_hat*T_I))*(1+(1/(2*CN0_L2_hat*T_I)));
            sigma2Measure_valueL5 = (1/(2*CN0_L5_hat*T_I))*(1+(1/(2*CN0_L5_hat*T_I)));
            Rk_adapt = diag([sigma2Measure_valueL1/2,sigma2Measure_valueL1/2,sigma2Measure_valueL2/2,sigma2Measure_valueL2/2,sigma2Measure_valueL5/2,sigma2Measure_valueL5/2]);
        end
    end
    
    Kk = Pkk1*JHk'*((JHk*Pkk1*JHk' + Rk_adapt)\eye(NumSeries)); % Kalman Gain computation
    xkk = xkk1 + Kk*(yk_formatted(k,:)' - ykk1);
    Pkk = Pkk1 - Kk*JHk*Pkk1;
    
%     gamma = 0.1;
%     if TopologySelector == 1
%         D = [zeros(1,L+M-1),zeros(1,NumSeries*P)];
%         d = ones(1,1)*gamma;
%         if xkk(L+M-1 + 1) <= 0
%             D(1,L+M-1 + 1) = 1;
%         else
%             D(1,L+M-1 + 1) = 0;
%         end
%     elseif TopologySelector == 2
%         D = [zeros(1,L+M-1),zeros(1,NumSeries*P);
%              zeros(1,L+M-1),zeros(1,NumSeries*P)];
%         d = ones(2,1)*gamma;
%         
%         if xkk(L+M-1 + 1) <= 0
%             D(1,L+M-1 + 1) = 1;
%         else
%             D(1,L+M-1 + 1) = 0;
%         end
%         
%         if xkk(L+M-1 + 2) <= 0
%             D(1,L+M-1 + 2) = 1;
%         else
%             D(1,L+M-1 + 2) = 0;
%         end
%         
%     elseif TopologySelector == 3
%         D = [zeros(1,L+M-1),zeros(1,NumSeries*P);
%              zeros(1,L+M-1),zeros(1,NumSeries*P);
%              zeros(1,L+M-1),zeros(1,NumSeries*P)];
%         d = ones(3,1)*gamma;
%         if xkk(L+M-1 + 1) <= 0
%             D(1,L+M-1 + 1) = 1;
%         else
%             D(1,L+M-1 + 1) = 0;
%         end
%         
%         if xkk(L+M-1 + 2) <= 0
%             D(1,L+M-1 + 2) = 1;
%         else
%             D(1,L+M-1 + 2) = 0;
%         end
%         
%         if xkk(L+M-1 + 3) <= 0
%             D(1,L+M-1 + 3) = 1;
%         else
%             D(1,L+M-1 + 3) = 0;
%         end
%     end
%     
%     if sum(sum(D,1),2) > 0
%         xkk = xkk - Pkk*D'*((D*Pkk*D' + (1e-9)*eye(TopologySelector))\eye(TopologySelector))*(D*xkk-d);
%     end
end

end
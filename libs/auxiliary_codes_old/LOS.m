% Project Title: EKF_VAR_Specific_Config

% Author: Rodrigo de Lima Florindo

% Version: 1.0

% Date: 22/06/2024

% Description: This code may be used to generate the state transition
% matrix and noise covariance matrix of the LOS dynamics for a generic
% multi-frequency system using the development proposed in section 2.1.

% Its output are:
% 1 - The state transition matrix related to the LOS dynamics processes -
% FDnew
% 2 - The noise covariance matrix related to the LOS dynamics - QDnew

function [FDnew,QDnew] = LOS(L,M,T,sigma,delta)

    syms T_I tau real
    syms s [M+L-1 1] real
    syms d [L,1] real

    Aux1 = sym(ones(L-1,1));

    for i=1:(L-1)
        Aux1(i,1) = 2*pi*d(i);
    end

    A1 = sym([zeros(L-1,1),Aux1,zeros(L-1,M-2)]);

    A2 = sym([[zeros(M-1,1),eye(M-1)];zeros(1,M)]);

    A2(1,2) = 2*pi*(d(L));

    Fw = sym([[zeros(L-1,L-1),A1];[zeros(M,L-1),A2]]);

    FD = sym(zeros(M+L-1,M+L-1));

    for j=0:(M-1)
        FD = FD + ((Fw*T_I)^(j))/factorial(j);
    end

    % Calculating the covariance QD

    Qe = sym(zeros(L+M-1,L+M-1,L+M-1));

    for i=1:(M+L-1)
        Qe(i,i,i) = s(i);
    end

    E = sym(zeros(M+L-1));

    for j=0:(M-1)
        E = E + ((Fw*(T_I-tau))^(j))/factorial(j);
    end

    QDi = sym(zeros(L+M-1,L+M-1,L+M-1));
    QD = sym(zeros(L+M-1,L+M-1));
    for i = 1:(M+L-1)
        expr = E * Qe(:,:,i) * E';
        QDi(:,:,i) = int(expr,tau,0,T_I);
        QD = QD + QDi(:,:,i);
    end

    %Substitution the values of T_I, s and d on FD and QD
    FDnew = FD; % Initialize FDnew
    for i = 1:length(d)
        FDnew = subs(FDnew, d(i), delta(i)); % Substitute each element of d with corresponding delta
    end

    FDnew = subs(FDnew, T_I, T);
    FDnew = double(FDnew);

    QDnew = QD; % Initialize FDnew
    for i = 1:length(d)
        QDnew = subs(QDnew, d(i), delta(i)); % Substitute each element of d with corresponding delta
    end

    for i = 1:length(sigma)
        QDnew = subs(QDnew, s(i), sigma(i)); % Substitute each element of d with corresponding delta
    end

    % Substitute T_I
    QDnew = subs(QDnew, T_I, T);
    QDnew = double(QDnew);
end
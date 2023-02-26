%% State-Space Model
% Dynamic Common Factor Model
% y(1t) = Lambda(1row)*B(t) + e1(t), e1(t) ~ iidN(0, sig2e1)
% y(2t) = Lambda(2row)*B(t) + e2(t), e2(t) ~ iidN(0, sig2e2)
% y(3t) = Lambda(3row)*B(t) + e3(t), e3(t) ~ iidN(0, sig2e3)
% y(4t) = Lambda(4row)*B(t) + e4(t), e4(t) ~ iidN(0, sig2e4)
% y(5t) = Lambda(5row)*B(t) + e5(t), e5(t) ~ iidN(0, sig2e5)
% y(6t) = Lambda(6row)*B(t) + e6(t), e6(t) ~ iidN(0, sig2e6)
% y(7t) = Lambda(7row)*B(t) + e7(t), e7(t) ~ iidN(0, sig2e7)
% y(8t) = Lambda(8row)*B(t) + e8(t), e8(t) ~ iidN(0, sig2e8)
% y(9t) = Lambda(9row)*B(t) + e9(t), e9(t) ~ iidN(0, sig2e9)
% y(10t) = Lambda(10row)*B(t) + e10(t), e10(t) ~ iidN(0, sig2e10)

% tau = [3 6 9 12 18 24 30 36 60 120]

% Measurement equation 
% Y(t) = L*B(t) + e(t) ~ iidN(0,R)

% Transition equation
% B(t) = [bL(t) bS(t) bC(t)]'
% B(t) = mu + F*B(t-1) + v(t), v(t) ~ iidN(0, Q) 


% SS Parameter
% C = [0 0 0 0 0 0 0 0 0 0]', 
% H = Lambda Matrix,
% R = diag([sig2e1; sig2e2; sig2e3; sig2e4; sig2e5; sig2e6; sig2e7; sig2e8;
% sig2e9; sig2e10])

% Mu = [Mu1 Mu2 Mu3]' 
% F = diag([phi1;phi2;phi3]) 
% Q = [sig2v1 Cov1 Cov2 ; Cov1 sig2v2 Cov3 ; Cov2 Cov3 sig2v3]

clear; 
clc; 
warning('off','all')
%% Step 1: DGP

data = readmatrix("Data.xlsx");
T = rows(data);     % Sample size 
k = 3;        % Three common factor

% initial Parameter

sig2e1 = 0.1; 
%sig2e2 = 0.2; 
%sig2e3 = 0.3;
%sig2e4 = 0.5;
%sig2e5 = 0.4;
%sig2e6 = 0.2;
%sig2e7 = 0.1;
%sig2e8 = 0.3;
%sig2e9 = 0.2;
%sig2e10 = 0.1;

mu1 = 0;
mu2 = 0;
mu3 = 0;

phi1 = 0.5;
phi2 = 0.3;
phi3 = 0.2;

sig2v1=0.4;
sig2v2=0.3;
sig2v3=0.2;

Cov1 = 0.1;
Cov2 = 0.1;
Cov3 = 0.1;

%tru_para = [log(sig2e1);log(sig2e2);log(sig2e3);log(sig2e4);log(sig2e5);log(sig2e6);log(sig2e7);log(sig2e8);log(sig2e9);log(sig2e10);mu1;mu2;mu3;phi1;phi2;phi3;log(sig2v1);log(sig2v2);log(sig2v3);log(Cov1);log(Cov2);log(Cov3)]; 



%% Step 2: Maxmimum Likelihood Estimation

% Block for each parameters

indR = 1;
indMu = [2;3;4];
indF = [5;6;7;];
indQ = [8;9;10;11;12;13];

% Structure variables
Sn.data = data;
Sn.indR = indR;
Sn.indMu = indMu;
Sn.indF = indF;
Sn.indQ = indQ;

% Initial values
%psi0 = [log(sig2e1);log(sig2e2);log(sig2e3);log(sig2e4);log(sig2e5);log(sig2e6);log(sig2e7);log(sig2e8);log(sig2e9);log(sig2e10);mu1;mu2;mu3;phi1;phi2;phi3;log(sig2v1);log(sig2v2);log(sig2v3);log(Cov1);log(Cov2);log(Cov3)];
psi0 = [log(sig2e1);mu1;mu2;mu3;phi1;phi2;phi3;log(sig2v1);log(sig2v2);log(sig2v3);log(Cov1);log(Cov2);log(Cov3)];

% Index
indbj = 1:rows(psi0);
indbj = indbj';
% printi = 1 => See the opimization produdure
% printi = 0 => NOT see the opimization produdure
printi = 1; 

% Optimization
[psimx, fmax,Vj, Vinv] = SA_Newton(@lnlik,@paramconst,psi0,Sn,printi,indbj);

% Estimates by Deltamethod
thetamx = maketheta(psimx,Sn);                  % Transform psi -> theta
grad = Gradpnew1(@maketheta,psimx,indbj,Sn);    % Gradient
cov_fnl = grad*Vj*grad';                        % Covariance Matrix
diag_cov = diag(cov_fnl);                       % Variance (diagonal)
stde = sqrt(diag_cov);                          % Standard deviation
t_val = thetamx./stde;                          % t-value
p_val = 2*(1 - cdf('t',abs(t_val),T-k));        % p-value

%% Step 4: Filtering and Smoothing
[C,H,R,Mu,F,Q] = makePara(thetamx,Sn);
[~, Beta_ttm, P_ttm]= KM_filter(C,H,R,Mu,F,Q,data);
Beta_LB = Beta_ttm - 1.95*sqrt(P_ttm); 
Beta_UB = Beta_ttm + 1.95*sqrt(P_ttm);

[Beta_tTm, P_tTm] = KM_smooth(C,H,R,Mu,F,Q,data);
Beta_LB_SM = Beta_tTm - 1.95*sqrt(P_tTm);
Beta_UB_SM = Beta_tTm + 1.95*sqrt(P_tTm);

%% Step 3: Table / Figure Results 
% Table
disp('===========================================================');
disp(['    Index ','  True Para ', ' Estimates ', ' t value ',  ' p value']);
disp('===========================================================');
disp([indbj tru_para thetamx t_val p_val]);
disp('===========================================================');

% Figure
% Index 
i = 1:rows(Beta_ttm); 

% Filtered value
figure
plot(i, cm ,'k', i, Beta_ttm, 'b:', 'LineWidth',1.5);
legend('True','Filtered'); 
title('True and Filtered Common Factor');

figure
plot(i, Beta_ttm ,'k', i, Beta_LB, 'b:', i, Beta_UB,'r:','LineWidth',1.5)
legend('Common Factor', 'Low Band', 'High Band');
title('Filtered Common Factor and Confidence Interval');

% Smoothed value
figure
plot(i, cm, 'k', i, Beta_tTm,'r:', 'LineWidth',1.5);
legend('True','Smoothed');
title('True and Smoothed Common Factor');

figure
plot(i, Beta_tTm ,'k', i, Beta_LB_SM, 'b:', i, Beta_UB_SM,'r:','LineWidth',1.5)
legend('Common Factor', 'Low Band', 'High Band');
title('Smoothed Common Factor and Confidence Interval');








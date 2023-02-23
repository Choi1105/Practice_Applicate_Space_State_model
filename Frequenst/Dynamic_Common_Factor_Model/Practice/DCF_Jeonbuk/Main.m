%% State-Space Model
% Dynamic Common Factor Model
% y(1t) = gam1*C(t) + e(1t) , e(1t) ~ iidN(0, sig2e1)
% y(2t) = gam2*C(t) + e(2t) , e(2t) ~ iidN(0, sig2e2)
% y(3t) = gam3*C(t) + e(3t) , e(3t) ~ iidN(0, sig2e3)
% y(4t) = gam4*C(t) + e(4t) , e(4t) ~ iidN(0, sig2e4)
% y(5t) = gam5*C(t) + e(5t) , e(5t) ~ iidN(0, sig2e5)
% y(6t) = gam6*C(t) + e(6t) , e(6t) ~ iidN(0, sig2e6)
% y(7t) = gam7*C(t) + e(7t) , e(7t) ~ iidN(0, sig2e7)

% C(t) = Mu + phi1*C(t-1) + phi2*C(t-2) + w(t), w(t) ~ iidN(0, sig2w) 

% Measurement equation
% y(t) = H*B(t) + w(t), w(t) ~ iidN(0,R) 

% Transition equation
% B(t) = Mu + F*B(t-1) + v(t), v(t) ~ iidN(0,Q) 

% SS Parameter
% C = [0 0 0 0 0 0 0]', 
% H = [gam1 0 ; gam2 0 ; gam3 0 ; gam4 0 ; gam5 0 ; gam6 0 ; gam7 0 ],
% R = [sig2e1 0 0 0 0 0 0 ; 0 sig2e2 0 0 0 0 0 ;  0 0 sig2e3 0 0 0 0 ;  0 0 0 sig2e4 0 0 0 ;  0 0 0 0 sig2e5 0 0 ; 0 0 0 0 0 sig2e6 0 ; 0 0 0 0 0 0 sig2e7 ;]

% Mu = [ mu 0 ]', 
% F = [ phi1 phi2 ; 1 0 ], 
% Q = [ sig2w 0 ; 0 0 ]

clear; 
clc; 
warning('off','all')
%% Step 1: DGP
% Data information

data = readmatrix("Jeonbuk_Data.csv");
T = rows(data);     % Sample size 
k = 1;        % One common factor

% initial Parameter
gam2 = 0.1;
gam3 = 0.1;
gam4 = 0.1;
gam5 = 0.1;
gam6 = 0.1;
gam7 = 0.1;

sig2e1 = 0.01; 
sig2e2 = 0.01; 
sig2e3 = 0.01;
sig2e4 = 0.01;
sig2e5 = 0.01;
sig2e6 = 0.01;
sig2e7 = 0.01;

% M.E.
mu = 0; 
phi1 = 0.01;
phi2 = 0.01;
sig2w=0.008;

tru_para = [gam2;gam3;gam4;gam5;gam6;gam7;sig2e1;sig2e2;sig2e3;sig2e4;sig2e5;sig2e6;sig2e7;mu;phi1;phi2;sig2w]; 


%% Step 2: Maxmimum Likelihood Estimation

% Block for each parameters
indH = [1;2;3;4;5;6];
indR = [7;8;9;10;11;12;13];
indMu = 14;
indF = [15;16];
indQ = 17;

% Structure variables
Sn.data = data;
Sn.indH = indH;
Sn.indR = indR;
Sn.indMu = indMu;
Sn.indF = indF;
Sn.indQ = indQ;

% Initial values
psi0 = [gam2;gam3;gam4;gam5;gam6;gam7;sig2e1;sig2e2;sig2e3;sig2e4;sig2e5;sig2e6;sig2e7;mu;phi1;phi2;sig2w]; 

% Index
indbj = 1:rows(psi0);
indbj = indbj';
% printi = 1 => See the opimization produdure
% printi = 0 => NOT see the opimization produdure
printi = 1; 

%% Optimization
[psimx, fmax,Vj, Vinv] = SA_Newton(@lnlik,@paramconst,psi0,Sn,printi,indbj); 

%% Estimates by Deltamethod
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
%figure
%plot(i, data ,'k', i, Beta_ttm, 'b:', 'LineWidth',1.5);
%legend('True','Filtered'); 
%title('True and Filtered Common Factor');

figure
plot(i, Beta_ttm ,'k', i, Beta_LB, 'b:', i, Beta_UB,'r:','LineWidth',1.5)
legend('Common Factor', 'Low Band', 'High Band');
title('Filtered Common Factor and Confidence Interval');

% Smoothed value
%figure
%plot(i, data, 'k', i, Beta_tTm,'r:', 'LineWidth',1.5);
%legend('True','Smoothed');
%title('True and Smoothed Common Factor');

figure
plot(i, Beta_tTm ,'k', i, Beta_LB_SM, 'b:', i, Beta_UB_SM,'r:','LineWidth',1.5)
legend('Common Factor', 'Low Band', 'High Band');
title('Smoothed Common Factor and Confidence Interval');








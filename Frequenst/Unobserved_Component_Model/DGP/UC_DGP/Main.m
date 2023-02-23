%% State-Space Model
% Unobserved Component Model
% Y(t) = X(t) + Z(t), 
% X(t) = mu + X(t-1) + u(t), u(t) ~ iidN (0, sig2u)
% Z(t) = phi1*Z(t-1) +phi2*Z(t-2) + v(t), a(t) ~ iidN(0, sig2v)

% Measurement equation
% y(t) = C + H*B(t) + e(t), e(t) ~ N(0,R) 

% Transition equation
% B(t) = Mu + F*B(t-1) + v(t), v(t) ~ N(0,Q) 

% SS Parameter
% C = 0, 
% H = ( 1  1  0 ), 
% R = 0

% Mu = (mu 0 0)', 
% F = [1 0 0 ; 0 phi1 phi2 ; 0 1 0], 
% Q = [sig2u  0 0;0 sig2v 0;0 0 0]

clear; 
clc; 

%% Step 1: Data Generating Process

% Data information
T = 1500;     % Sample size 
cut = 500;    % Burn-in size

% True Parameter
% AR
phi1 = 0.5;
phi2 = 0.2;
sig2v = 0.5;

% RW
mu = 0.2;
sig2u = 0.1;

% Set of true parameters
tru_para = [phi1;phi2;sig2u;sig2v;mu];

% Pre-allocation for x(t), z(t)
xm = zeros(T,1); 
zm = zeros(T,1);

% Initial x(t), z(t) for t = 1, 2
xm(1) = mu + randn(1,1)*sqrt(sig2u);
xm(2) = mu + xm(1) + randn(1,1)*sqrt(sig2u);

zm(1) = randn(1,1)*sqrt(sig2v);
zm(2) = phi1*zm(1)+randn(1,1)*sqrt(sig2v);

% Generate x(t), z(t) for t = 3 ~ T
for t = 3:T
    
    xm(t) = mu + xm(t-1) + randn(1,1)*sqrt(sig2u);
    zm(t) = phi1*zm(t-1) + phi2*zm(t-2) + randn(1,1)*sqrt(sig2v);
    
end

% Generate y(t)
ym = xm + zm;

% Burn-in
ym = ym(cut+1:end);
xm = xm(cut+1:end);
zm = zm(cut+1:end);

%% Step 2: Maxmimum Likelihood Estimation
% Data
data = ym;

% Block for each parameters
indF = [1;2];
indQ = [3;4];
indMu = 5;

% Structure variables
Sn.data = data;
Sn.indF = indF;
Sn.indQ = indQ;
Sn.indMu = indMu;

% Initial values
psi0 = [phi1;phi2;log(sig2u);log(sig2v);mu];

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
p_val = 2*(1 - cdf('t',abs(t_val),T-5));        % p-value

%% Step 4: Filtering and Smoothing
[C,H,R,Mu,F,Q] = makePara(thetamx,Sn);
[~, Beta_ttm, P_ttm]= KM_filter(C,H,R,Mu,F,Q,data);
Beta_LB = Beta_ttm - 1.95*sqrt(P_ttm); 
Beta_UB = Beta_ttm + 1.95*sqrt(P_ttm);

[Beta_tTm, P_tTm] = KM_smooth(C,H,R,Mu,F,Q,data);
Beta_LB_SM = Beta_tTm - 1.95*sqrt(P_tTm);
Beta_UB_SM = Beta_tTm + 1.95*sqrt(P_tTm);

% Burn some samples b/c wild guess
burn = 3;
ym = ym(burn+1:end);
xm = xm(burn+1:end);
zm = zm(burn+1:end);
Beta_ttm = Beta_ttm(burn+1:end,:);
P_ttm = P_ttm(burn+1:end,:); 
Beta_LB = Beta_LB(burn+1:end,:);
Beta_UB = Beta_UB(burn+1:end,:);
Beta_tTm = Beta_tTm(burn+1:end,:);
P_tTm = P_tTm(burn+1:end,:);
Beta_LB_SM = Beta_LB_SM(burn+1:end,:);
Beta_UB_SM = Beta_UB_SM(burn+1:end,:);

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

% Filtered values
figure
plot(i, Beta_ttm(:,1) ,'k', i, Beta_LB(:,1), 'b:', i, Beta_UB(:,1),'r:','LineWidth',1.5)
legend('Trend', 'Low Band', 'High Band');
title('Filtered Trend and Confidence Interval');

figure
plot(i, Beta_ttm(:,2) ,'k', i, Beta_LB(:,2), 'b:', i, Beta_UB(:,2),'r:','LineWidth',1.5)
legend('Cycle', 'Low Band', 'High Band');
title('Filtered Cycle and Confidence Interval');

figure
plot(i, xm ,'k', i, Beta_ttm(:,1), 'b:', 'LineWidth',1.5);
legend('True','Filtered'); 
title('True and Filtered Trend');

figure
plot(i, zm ,'k', i, Beta_ttm(:,2), 'r:', 'LineWidth',1.5);
legend('True','Filtered'); 
title('True and Filtered Cycle');

figure
plot(i, ym ,'k', i, Beta_ttm(:,1), 'b:', i, Beta_ttm(:,2), 'r:', 'LineWidth',1.5);
legend('Data','Trend','Cycle'); 
title('Data and Filtered Trend, Cycle');

% Smoothed values
figure
plot(i, Beta_tTm(:,1) ,'k', i, Beta_LB_SM(:,1), 'b:', i, Beta_UB_SM(:,1),'r:','LineWidth',1.5)
legend('Trend', 'Low Band', 'High Band');
title('Smoothed Trend and Confidence Interval');

figure
plot(i, Beta_tTm(:,2) ,'k', i, Beta_LB_SM(:,2), 'b:', i, Beta_UB_SM(:,2),'r:','LineWidth',1.5)
legend('Cycle', 'Low Band', 'High Band');
title('Smoothed Cycle and Confidence Interval');

figure
plot(i, xm, 'k', i, Beta_tTm(:,1) ,'b:', 'LineWidth',1.5);
legend('True','Smoothed');
title('True and Smoothed Trend');

figure
plot(i, zm, 'k', i, Beta_tTm(:,2) ,'r:', 'LineWidth',1.5);
legend('True','Smoothed');
title('True and Smoothed Cycle');

figure
plot(i, ym ,'k', i, Beta_tTm(:,1), 'b:', i, Beta_tTm(:,2),'r:' , 'LineWidth',1.5);
legend('Data','Trend','Cycle'); 
title('Data and Smoothed Trend, Cycle');
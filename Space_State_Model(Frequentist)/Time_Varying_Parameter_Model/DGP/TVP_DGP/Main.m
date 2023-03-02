%% State-Space Model
% Time Varying Parameter Model

% Y(t) = X(t)*B(t) + e(t), e(t) ~ iidN (0, sig2e) 
% B(t) = phi*B(t-1) + v(t), v(t) ~ iidN(0, sig2v)
% X(t) exogenous variable

% Measurement equation
% y(t) = C + H*B(t) + e(t), e(t) ~ N(0,R) 

% Transition equation
% B(t) = Mu + F*B(t-1) + v(t), v(t) ~ N(0,Q) 

% SS Parameter
% C = 0, 
% H = X(t), 
% R = sig2e

% Mu = 0, 
% F = 1, 
% Q = sig2v

clear;
clc; 

%% Step 1: Data Generating Process

% Data information
T = 1500;     % Sample size 
cut = 500;    % Burn-in size
k = 1;        % # of latent variables
              % (Actual sample size = T - cut)
              
% True Parameter
sig2e = 0.5;
sig2v = 1; 

% Set of true parameters
tru_para = [sig2e;sig2v];

% Given exogenous variable
xm = 5*randn(T,1);  

% Pre-allocation for y(t), b(t)
ym = zeros(T,1); 
bm = zeros(T,1);

% Initial y(t), b(t) for t = 1
bm(1) = randn(1,1)*sqrt(sig2v);
ym(1) = xm(1)*bm(1) + randn(1,1)*sqrt(sig2e);

for t = 2:T
    
    bm(t) = bm(t-1) + randn(1,1)*sqrt(sig2v);
    ym(t) = xm(t)*bm(t) + randn(1,1)*sqrt(sig2e);
    
end

% Burn-in
ym = ym(cut+1:end);
xm = xm(cut+1:end);
bm = bm(cut+1:end);


%% Step 2: Maxmimum Likelihood Estimation
% Block for each parameters
indR = 1;
indQ = 2;

% Structure variables
Sn.ym = ym;
Sn.xm = xm;
Sn.indR = indR;
Sn.indQ = indQ;

% Initial values
psi0 = [log(sig2e);log(sig2v)];

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
[~, Beta_ttm, P_ttm]= KM_filter(C,H,R,Mu,F,Q,ym);
Beta_LB = Beta_ttm - 1.95*sqrt(P_ttm); 
Beta_UB = Beta_ttm + 1.95*sqrt(P_ttm);

[Beta_tTm, P_tTm] = KM_smooth(C,H,R,Mu,F,Q,ym);
Beta_LB_SM = Beta_tTm - 1.95*sqrt(P_tTm);
Beta_UB_SM = Beta_tTm + 1.95*sqrt(P_tTm);

% Burn some samples b/c wild guess
burn = 3;
ym = ym(burn+1:end);
bm = bm(burn+1:end);
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
plot(i, Beta_ttm, 'k', i, Beta_LB, 'b:', i, Beta_UB,'r:','LineWidth',1.5)
legend('Time-varying Parameters', 'Low Band', 'High Band');
title('Filtered Time-varying Parameters');


figure
plot(i, bm ,'k', i, Beta_ttm, 'b:', 'LineWidth',1.5);
legend('Time-varying Parameters','Filtered'); 
title('True and Filtered Time-varying Parameters');

% Smoothed values
figure
plot(i, Beta_tTm ,'k', i, Beta_LB_SM, 'b:', i, Beta_UB_SM,'r:','LineWidth',1.5)
legend('Time-varying Parameters', 'Low Band', 'High Band');
title('Smoothed Time-varying Parameters and Confidence Interval');


figure
plot(i, bm, 'k', i, Beta_tTm,'b:', 'LineWidth',1.5);
legend('True','Smoothed');
title('True and Smoothed Time-varying Parameters');


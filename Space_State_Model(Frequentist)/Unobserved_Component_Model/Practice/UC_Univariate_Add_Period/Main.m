

%% State-Space Model

% Unobserved Component Model (Bivariate)
% y(t) = n(t) + x(t),
% n(t) = Mu + n(t-1) + v(t), v(t) ~ iidN (0, sig2v)
% x(t) = phi1*x(t-1) +phi2*x(t-2) + e(t), e(t) ~ iidN(0, sig2e)

% Where y(t) = log of Real GDP (US Real Gross Domestic Product - FRED)
%       n(t) = Stochastic trend component 
%       x(t) = Stationary cyclical component
%       v(t), e(t) = independent white noise processes 

% 1. Measurement equation
%    y(t) = H*B(t)

% 1. Transition equation
%    B(t) = Mu + F*B(t-1) + u(t), u(t) ~ N(0,Q) 

% 1. SS Parameter
%    C = 0, 
%    H = ( 1  1  0 ), 
%    R = 0 

% Mu = [ Mu ; 0 ; 0 ]
% F = [ 1 0 0 ; 0 phi1 phi2 ; 0 1 0 ], 
% Q = [ sig2v 0 0 ; 0 sig2e 0 ; 0 0 0 ]

clear; 
clc; 
%% Data
data = log(readmatrix("Real_GDP_Add_Period.xlsx", "Range","A1:A175"));
ym = data;
T = rows(data);

% Initial Parameter
% AR
phi1 = 0.7;
phi2 = -0.2;
sig2e = 0.01;
sig2v = 0.01;
mu = 0;

%% Step 2: Maxmimum Likelihood Estimation

% Block for each parameters
indF = [1;2];
indQ = [3;4];
indMu = 5;

% Structure variables
Sn.data = data;
Sn.indF = indF;
Sn.indQ = indQ;
Sn.indMu = indMu;

psi0 = [phi1;phi2;sig2v;sig2e;mu];

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
data = data(burn+1:end);
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
disp([indbj psi0 thetamx t_val p_val]);
disp('===========================================================');

% Figure
% Index 
i = 1:rows(Beta_ttm); 

% Filtered values
tiledlayout(2,3)
nexttile
plot(i, Beta_ttm(:,1) ,'k')%, i, Beta_LB(:,1), 'b:', i, Beta_UB(:,1),'r:','LineWidth',1.5)
legend('Trend')%, 'Low Band', 'High Band');
title('Filtered Trend and Confidence Interval');

nexttile
plot(i, Beta_ttm(:,2) ,'k')%, i, Beta_LB(:,2), 'b:', i, Beta_UB(:,2),'r:','LineWidth',1.5)
legend('Cycle')%, 'Low Band', 'High Band');
title('Filtered Cycle')% and Confidence Interval');

nexttile
plot(i, data ,'k', i, Beta_ttm(:,1), 'b:', i, Beta_ttm(:,2), 'r:', 'LineWidth',1.5);
legend('Data','Trend','Cycle'); 
title('Data and Filtered Trend, Cycle');

% Smoothed values
nexttile
plot(i, Beta_tTm(:,1) ,'k')%, i, Beta_LB_SM(:,1), 'b:', i, Beta_UB_SM(:,1),'r:','LineWidth',1.5)
legend('Trend')%, 'Low Band', 'High Band');
title('Smoothed Trend')% and Confidence Interval');

nexttile
plot(i, Beta_tTm(:,2) ,'k')%, i, Beta_LB_SM(:,2), 'b:', i, Beta_UB_SM(:,2),'r:','LineWidth',1.5)
legend('Cycle')%, 'Low Band', 'High Band');
title('Smoothed Cycle')% and Confidence Interval');

nexttile
plot(i, data ,'k', i, Beta_tTm(:,1), 'b:', i, Beta_tTm(:,2),'r:' , 'LineWidth',1.5);
legend('Data','Trend','Cycle'); 
title('Data and Smoothed Trend, Cycle');
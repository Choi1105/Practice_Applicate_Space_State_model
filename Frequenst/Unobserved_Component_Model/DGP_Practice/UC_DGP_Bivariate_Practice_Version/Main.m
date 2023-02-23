%% State-Space Model
% This code sheet : For applicate Real Data

% Unobserved Component Model (Bivariate)
% y(t) = n(t) + x(t),
% n(t) = g(t-1) + n(t-1) + v(t), v(t) ~ iidN (0, sig2v)
% g(t) = g(t-1) + w(t), w(t) ~ iidN (0, sig2w)
% x(t) = phi1*x(t-1) +phi2*x(t-2) + e(t), e(t) ~ iidN(0, sig2e)
% U(t) = L(t) + C(t)
% L(t) = L(t-1) + vl(t), vl(t) ~ iidN (0, sig2vl)
% C(t) = a0*x(t) + a1*x(t-1) + a2*x(t-2) + ec(t), ec(t) ~ iidN(0, sig2ec)
 

% Where y(t) = log of Real GDP
%       n(t) = Stochastic trend component 
%       x(t) = Stationary cyclical component
%       v(t), w(t), e(t) = independent white noise processes 
%       g(t) = Stochastic trend component's drift term
%       L(t) = trend component
%       C(t) = stationary component(실질생산의 cyclical component(x(t))
%       function)

% 1. Measurement equation
%    y(t)  = H*B(t) + ey(t)

% 1. Transition equation
%    B(t) = F*B(t-1) + u(t), u(t) ~ N(0,Q) 

% 1. SS Parameter
%    C = 0, 
%    H = [ 1 0 0 1 0 0 ; 0 0 1 a0 a1 a2 ], 
%    R = [ 0 ; ec(t) ] 

%    Mu = [ 0 0 0 0 0 0 ]', 
%    F = [1 1 0 0 0 0 ; 0 1 0 0 0 0 ;  0 0 1 0 0 0  ; 0 0 0 phi1 phi2 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0], 
%    Q = [sig2v 0 0 0 0 0; 0 sig2w 0 0 0 0 ; 0 0 sig2vl 0 0 0 ; 0 0 0 sig2e 0 0 ; 0 0 0 0 0 0 ; 0 0 0 0 0 0 ]

clear; 
clc; 

%% Step 1: Data Generating Process
T = 1500;     % Sample size

% True Parameter
% AR
phi1 = 0.5;
phi2 = 0.3;
sig2e = 0.1;

% RW
mu = 0;
sig2w = 0.2;

% Stochastic Trend component
sig2v = 0.2;
sig2vl = 0.1;
sig2ec = 0.3;

alp0 = 0.4;
alp1 = 0.2;
alp2 = 0.1;

% Set of true parameters
tru_para = [phi1;phi2;sig2v;sig2w;sig2vl;sig2e;alp0;alp1;alp2;sig2ec;mu];

% Pre-allocation for n(t), g(t), x(t)
nm = zeros(T,1); 
gm = zeros(T,1);
xm = zeros(T,1);
Lm = zeros(T,1);
Cm = zeros(T,1);

% Initial xm(t), zm(t) for t = 1, 2
gm(1) = randn(1,1)*sqrt(sig2w);
gm(2) = gm(1) + randn(1,1)*sqrt(sig2w);

nm(1) = randn(1,1)*sqrt(sig2w) + randn(1,1)*sqrt(sig2v);
nm(2) = gm(1) + nm(1) + randn(1,1)*sqrt(sig2v);

xm(1) = randn(1,1)*sqrt(sig2e);
xm(2) = phi1*xm(1) + randn(1,1)*sqrt(sig2e);
xm(3) = phi1*xm(2) + phi2*xm(1) + randn(1,1)*sqrt(sig2e);

Lm(1) = randn(1,1)*sqrt(sig2vl);
Lm(2) = Lm(1) + randn(1,1)*sqrt(sig2vl);

Cm(1) = alp0*xm(1) + randn(1,1)*sqrt(sig2ec);
Cm(2) = alp0*xm(2) + alp1*xm(1) + randn(1,1)*sqrt(sig2ec);
Cm(3) = alp0*xm(3) + alp1*xm(2) + alp2*xm(1) + randn(1,1)*sqrt(sig2ec);

% Generate nm(t), gm(t), xm(t) for t = 3 ~ T
for t = 3:T
    gm(t) = gm(t-1) + randn(1,1)*sqrt(sig2w);
    nm(t) = gm(t-1) + nm(t-1) + randn(1,1)*sqrt(sig2v);
    xm(t) = phi1*xm(t-1) + phi2*xm(t-2) + randn(1,1)*sqrt(sig2e);
    Lm(t) = Lm(t-1) + randn(1,1)*sqrt(sig2vl);
end

for t = 4:T
    Cm(t) = alp0*xm(t) + alp1*xm(t-1) + alp2*xm(t-2) +randn(1,1)*sqrt(sig2ec);
end

% Generate 
Ym = nm + xm;
Um = Lm + Cm;
ym = [Ym Um];

%% Step 2: Maxmimum Likelihood Estimation
% Data
data = ym;

% Block for each parameters
indF = [1;2];
indQ = [3;4;5;6];
indH = [7;8;9];
indR = 10;
indMu = 11;

% Structure variables
Sn.data = data;
Sn.indF = indF;
Sn.indQ = indQ;
Sn.indMu = indMu;
Sn.indR = indR;
Sn.indH = indH;

% Initial values
psi0 = [phi1;phi2;log(sig2v);log(sig2w);log(sig2vl);log(sig2e);alp0;alp1;alp2;sig2ec;mu];

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
plot(i, Lm ,'k', i, Beta_ttm(:,1), 'b:', 'LineWidth',1.5);
legend('True','Filtered'); 
title('True and Filtered Trend');

figure
plot(i, Cm ,'k', i, Beta_ttm(:,2), 'r:', 'LineWidth',1.5);
legend('True','Filtered'); 
title('True and Filtered Cycle');

figure
plot(i, Um ,'k', i, Beta_ttm(:,1), 'b:', i, Beta_ttm(:,2), 'r:', 'LineWidth',1.5);
legend('Data','Trend','Cycle'); 
title('Data and Filtered Trend, Cycle');

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
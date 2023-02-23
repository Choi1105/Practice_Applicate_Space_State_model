%% State-Space Model
% Unobserved Component Model(bivariate)
% Y(t) = [y(t), U(t)]'
% n(t) = g(t-1) + n(t-1) + v(t), v(t) ~ iidN (0, sig2v)
% g(t) = g(t-1) + w(t), w(t) ~ iidN(0, sig2w)
% l(t) = l(t-1) + vl(t), vl(t) ~ iidN(0, sig2vl)
% x(t) = phi1*x(t-1) +phi2*x(t-2) + e(t), e(t) ~ iidN(0, sig2e)
% c(t) = a0*x(t) + a1*x(t-1) + a2*x(t-2) + c(t) ~ iidN(0, sig2c)

% Measurement equation
% Y(t) = C + H*B(t) + a(t), (t) ~ N(0,R) 

% Transition equation
% B(t) = Mu + F*B(t-1) + u(t), u(t) ~ N(0,Q) 

% SS Parameter
% C = 0, 
% H = ( 1 0 0 1 0 0 ; 0 0 1 a0 a1 a2 ), 
% R = ( 0 0 ; 0 sig2c )

% Mu = (0 0 0 0 0 0)', 
% F = [1 1 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 phi1 phi2 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0], 
% Q = [sig2v 0 0 0 0 0 ; 0 sig2w 0 0 0 0 ; 0 0 sig2vl 0 0 0 ; 0 0 0 sig2e 0 0 ; 0 0 0 0 0 0 ; 0 0 0 0 0 0]

clear; 
clc; 

%% Step 1: Data Generating Process

% Data information
T = 1500;     % Sample size 
cut = 500;    % Burn-in size

% True Parameter (Real GDP)
% AR
phi1 = 0.5;
phi2 = 0.2;
sig2e = 0.5;

% RW
mu = 0;
sig2w = 0.1;

% Permanent Component(Trend)
sig2v = 0.2;

% True Parameter (Unemployment Rate)
% RW
sig2vl = 0.2;

% AR
a0 = 0.4;
a1 = 0.2;
a2 = 0.1;
sig2c = 0.2;

% Set of true parameters
tru_para = [phi1;phi2;sig2v;sig2w;sig2vl;sig2e;a0;a1;a2;sig2c;mu];

% Pre-allocation for N(t), G(t), X(t)
nm = zeros(T,1);
gm = zeros(T,1);
xm = zeros(T,1);
lm = zeros(T,1);
cm = zeros(T,1);

% Initial n(t), g(t), x(t), l(t), c(t) for t = 1, 2
gm(1) = randn(1,1)*sqrt(sig2w);
gm(2) = gm(1) + randn(1,1)*sqrt(sig2w);

nm(1) = randn(1,1)*sqrt(sig2v);
nm(2) = gm(1) + nm(1) + randn(1,1)*sqrt(sig2v);

xm(1) = randn(1,1)*sqrt(sig2e);
xm(2) = phi1*xm(1)+randn(1,1)*sqrt(sig2e);

lm(1) = randn(1,1)*sqrt(sig2vl);
lm(2) = lm(1) + randn(1,1)*sqrt(sig2vl);

cm(1) = a0*xm(1) + randn(1,1)*sqrt(sig2c);
cm(2) = a0*xm(2) + a1*xm(1) + randn(1,1)*sqrt(sig2c);

% Generate Nm(t), Gm(t), Xm(t) for t = 3 ~ T
for t = 3:T
    
    nm(t) = gm(t-1) + nm(t-1) + randn(1,1)*sqrt(sig2v);
    gm(t) = gm(t-1) + randn(1,1)*sqrt(sig2w);
    xm(t) = phi1*xm(t-1) + phi2*xm(t-2) + randn(1,1)*sqrt(sig2e);
    lm(t) = lm(t-1) + randn(1,1)*sqrt(sig2vl);
    cm(t) = a0*xm(t) + a1*xm(t-1) + a2*xm(t-2) + randn(1,1)*sqrt(sig2c);
    
end

% Generate Y(t)
ym = nm + xm;
Um = lm + cm;

% Burn-in
ym = ym(cut+1:end);
nm = nm(cut+1:end);
gm = gm(cut+1:end);
xm = xm(cut+1:end);
Um = Um(cut+1:end);
lm = lm(cut+1:end);
cm = cm(cut+1:end);
Ym = [ym, Um];

%% Step 2: Maxmimum Likelihood Estimation
% Data
data = Ym;

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
Sn.indH = indH;
Sn.indR = indR;
Sn.indMu = indMu;

% Initial values
psi0 = [phi1;phi2;log(sig2v);log(sig2w);log(sig2vl);log(sig2e);a0;a1;a2;log(sig2c);mu];

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
burn = 5;
ym = ym(burn+1:end);
nm = nm(burn+1:end);
gm = gm(burn+1:end);
xm = xm(burn+1:end);
Um = Um(burn+1:end);
lm = lm(burn+1:end);
cm = cm(burn+1:end);
Ym = [ym, Um];

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
plot(i, nm ,'k', i, Beta_ttm(:,1), 'b:', 'LineWidth',1.5);
legend('True','Filtered'); 
title('True and Filtered Trend');

figure
plot(i, xm ,'k', i, Beta_ttm(:,2), 'r:', 'LineWidth',1.5);
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
plot(i, nm, 'k', i, Beta_tTm(:,1) ,'b:', 'LineWidth',1.5);
legend('True','Smoothed');
title('True and Smoothed Trend');

figure
plot(i, xm, 'k', i, Beta_tTm(:,2) ,'r:', 'LineWidth',1.5);
legend('True','Smoothed');
title('True and Smoothed Cycle');

figure
plot(i, ym ,'k', i, Beta_tTm(:,1), 'b:', i, Beta_tTm(:,2),'r:' , 'LineWidth',1.5);
legend('Data','Trend','Cycle'); 
title('Data and Smoothed Trend, Cycle');
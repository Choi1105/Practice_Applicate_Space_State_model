%% State-Space Model
% Time Varying Parameter Model

% Y(t) = B0(t) + B1(t)*Pi(t-1) + B2(t)*g(t-1) + B3(t)*i(t-1) + e(t), e(t) ~ iidN (0, sig2e) 
% B0(t) = B0(t-1) + v0(t), v0(t) ~ iidN(0, sig2v0)
% B1(t) = B1(t-1) + v1(t), v1(t) ~ iidN(0, sig2v1)
% B2(t) = B2(t-1) + v2(t), v2(t) ~ iidN(0, sig2v2)
% B3(t) = B3(t-1) + v3(t), v3(t) ~ iidN(0, sig2v3)

% Pi(t-1), g(t-1), i(t-1) are exogenous variable

% Measurement equation
% y(t) = H*B(t) + e(t), e(t) ~ N(0,R) 

% Transition equation
% Bi(t) = Bi(t-1) + vi(t), vi(t) ~ N(0,Q) i = 0,1,2,3
 

% SS Parameter
% C = 0, 
% H = [1 pi(t-1) g(t-1) i(t-1)]
% R = sig2e

% Mu = [0 0 0 0]' 
% F = eye(4) 
% Q = [sig2v0 0 0 0 ; 0 sig2v1 0 0 ; 0 0 sig2v2 0 ; 0 0 0 sig2v3]

clear;
clc; 

%% Step 1: Data Generating Process
% Data information
T = 2000;     % Sample size 
cut = 500;    % Burn-in size
k = 1;        % # of latent variables
              % (Actual sample size = T - cut)
              
% True Parameter
sig2e = 0.14;

sig2v0 = 0.22; 
sig2v1 = 0.13; 
sig2v2 = 0.14; 
sig2v3 = 0.17; 

% Set of true parameters
tru_para = [sig2e;sig2v0;sig2v1;sig2v2;sig2v3];

% Given exogenous variable
xm0 = ones(T,1);
xm1 = 5*randn(T,1); 
xm2 = 5*randn(T,1);  
xm3 = 5*randn(T,1);  

% Pre-allocation for y(t), b(t)
ym = zeros(T,1); 
bm0 = zeros(T,1);
bm1 = zeros(T,1);
bm2 = zeros(T,1);
bm3 = zeros(T,1);

bm = [bm0 bm1 bm2 bm3]';
xm = [xm0 xm1 xm2 xm3];

% Initial y(t), b(t) for t = 1
bm0(1) = randn(1,1)*sqrt(sig2v0);
bm1(1) = randn(1,1)*sqrt(sig2v1);
bm2(1) = randn(1,1)*sqrt(sig2v2);
bm3(1) = randn(1,1)*sqrt(sig2v3);

ym(1) = xm(1)*bm(1) + randn(1,1)*sqrt(sig2e);

for t = 2:T
    
    bm0(t) = bm0(t-1) + randn(1,1)*sqrt(sig2v0);
    bm1(t) = bm1(t-1) + randn(1,1)*sqrt(sig2v1);
    bm2(t) = bm2(t-1) + randn(1,1)*sqrt(sig2v2);
    bm3(t) = bm3(t-1) + randn(1,1)*sqrt(sig2v3);

    ym(t) = xm0(t-1)*bm0(t) + xm1(t-1)*bm1(t) + xm2(t-1)*bm2(t) + xm3(t-1)*bm3(t) + randn(1,1)*sqrt(sig2e);
    
end


% Burn-in
ym = ym(cut+1:end);
xm0 = xm0(cut+1:end);
bm0 = bm0(cut+1:end);
xm1 = xm1(cut+1:end);
bm1 = bm1(cut+1:end);
xm2 = xm2(cut+1:end);
bm2 = bm2(cut+1:end);
xm3 = xm3(cut+1:end);
bm3 = bm3(cut+1:end);

bm = [bm0 bm1 bm2 bm3];
xm = [xm0 xm1 xm2 xm3];

%% Step 2: Maxmimum Likelihood Estimation
% Block for each parameters
indR = 1;
indQ = [2;3;4;5];

% Structure variables
Sn.ym = ym;
Sn.xm = xm;
Sn.indR = indR;
Sn.indQ = indQ;

% Initial values
psi0 = [sig2e;sig2v0;sig2v1;sig2v2;sig2v3];

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
tiledlayout(4,4)
for k = 1:4
    nexttile
    plot(i, Beta_ttm(:,k), 'k', i, Beta_LB(:,k), 'b:', i, Beta_UB(:,k),'r:','LineWidth',1.5)
    legend('Time-varying Parameters', 'Low Band', 'High Band');
    title('Filtered Time-varying Parameters');
    nexttile
    plot(i, bm(:,k) ,'k', i, Beta_ttm(:,k), 'b:', 'LineWidth',1.5); 
    legend('Time-varying Parameters','Filtered'); 
    title('True and Filtered Time-varying Parameters');
end

% Smoothed values

for j = 1:4
    nexttile
    plot(i, Beta_tTm(:,j) ,'k', i, Beta_LB_SM(:,j), 'b:', i, Beta_UB_SM(:,j),'r:','LineWidth',1.5)
    legend('Time-varying Parameters', 'Low Band', 'High Band');
    title('Smoothed Time-varying Parameters and Confidence Interval');

    nexttile
    plot(i, bm(:,j), 'k', i, Beta_tTm(:,j),'b:', 'LineWidth',1.5);
    legend('True','Smoothed');
    title('True and Smoothed Time-varying Parameters');
end

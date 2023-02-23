%% State-Space Model
% Time Varying Parameter Model

% y(t) = X(t-1)*B(t) + e(t), e(t) ~ iidN (0, sig2e) 
% Bi(t) = Bi(t-1) + vi(t), vi(t) ~ iidN(0, sig2vi) i = 0, 1, 2, 3
% X(t-1) exogenous variable
% X(t-1) = [ 1 pi(t-1) g(t-1) i(t-1) ]
% Bi(t) = [ b0(t) b1(t) b2(t) b3(t) ]' 

% Measurement equation
% y(t) = X(t-1)*B(t) + e(t), e(t) ~ N(0,R) 

% Transition equation
% B(t) = B(t-1) + v(t), v(t) ~ N(0,Q) 
 

% SS Parameter
% C = 0, 
% H = X(t-1),
% R = sig2e

% Mu = [0 0 0 0]', 
% F = eye(4), 
% Q = [sig2v1 0 0 0 ; 0 sig2v2 0 0 ; 0 0 sig2v3 0 ; 0 0 0 sig2v4]
clear;
clc; 

%% Step 1: Data Generating Process

% Data information
ym = readmatrix("TVP_DATA.xlsx", 'Range', 'C2:C261');
x2m = readmatrix("TVP_DATA.xlsx", 'Range', 'A1:A260');
x3m = readmatrix("TVP_DATA.xlsx", 'Range', 'B1:B260');
x4m = readmatrix("TVP_DATA.xlsx", 'Range', 'C1:C260');
T = rows(ym);
x1m = ones(T,1);
k = 1;

xm = [x1m x2m x3m x4m];


% Initial Parameter
sig2e = 0.2123;
sig2v0 = 0.06322;
sig2v1 = 0.0345;
sig2v2 = 0.08234;
sig2v3 = 0.0567875;


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
figure
plot(i, Beta_ttm, 'k', i, Beta_LB, 'b:', i, Beta_UB,'r:','LineWidth',1.5)
legend('Time-varying Parameters', 'Low Band', 'High Band');
title('Filtered Time-varying Parameters');


%figure
%plot(i, bm ,'k', i, Beta_ttm, 'b:', 'LineWidth',1.5);
%legend('Time-varying Parameters','Filtered'); 
%title('True and Filtered Time-varying Parameters');

% Smoothed values
figure
plot(i, Beta_tTm ,'k', i, Beta_LB_SM, 'b:', i, Beta_UB_SM,'r:','LineWidth',1.5)
legend('Time-varying Parameters', 'Low Band', 'High Band');
title('Smoothed Time-varying Parameters and Confidence Interval');


%figure
%plot(i, bm, 'k', i, Beta_tTm,'b:', 'LineWidth',1.5);
%legend('True','Smoothed');
%title('True and Smoothed Time-varying Parameters');


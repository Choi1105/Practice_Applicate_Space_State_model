%% State-Space Model
% Time Varying Parameter Model

% Y(t) = B1(t)*Pi(t-1) + B2(t)*g(t-1) + B3(t)*i(t-1) + e(t), e(t) ~ iidN (0, sig2e) 
% B1(t) = B1(t-1) + v1(t), v1(t) ~ iidN(0, sig2v1)
% B2(t) = B2(t-1) + v2(t), v2(t) ~ iidN(0, sig2v2)
% B3(t) = B3(t-1) + v3(t), v3(t) ~ iidN(0, sig2v3)

% Pi(t-1), g(t-1), i(t-1) are exogenous variable

% Measurement equation
% y(t) = H*B(t) + e(t), e(t) ~ N(0,R) 

% Transition equation
% Bi(t) = Bi(t-1) + vi(t), vi(t) ~ N(0,Q) i = 1,2,3
 

% SS Parameter
% C = 0, 
% H = [pi(t-1) g(t-1) i(t-1)]
% R = sig2e

% Mu = [0 0 0 ]' 
% F = eye(3) 
% Q = [sig2v1 0 0 ; 0 sig2v2 0 ; 0 0 sig2v3]

clear;
clc; 

%% Step 1: Data Generating Process

% Data information
ym = readmatrix("TVP_DATA.xlsx", 'Range', 'C2:C261');
xm1 = readmatrix("TVP_DATA.xlsx", 'Range', 'A1:A260');
xm2 = readmatrix("TVP_DATA.xlsx", 'Range', 'B1:B260');
xm3 = readmatrix("TVP_DATA.xlsx", 'Range', 'C1:C260');
T = rows(ym);
k = 1;

xm = [xm1 xm2 xm3];


% Initial Parameter
sig2e = 0.2123;

sig2v1 = 0.0345;
sig2v2 = 0.08234;
sig2v3 = 0.0567875;


%% Step 2: Maxmimum Likelihood Estimation
% Block for each parameters
indR = 1;
indQ = [2;3;4];

% Structure variables
Sn.ym = ym;
Sn.xm = xm;
Sn.indR = indR;
Sn.indQ = indQ;

% Initial values
psi0 = [sig2e;sig2v1;sig2v2;sig2v3];

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
tiledlayout(3,2);
for m = 1:3
    nexttile
    plot(i, Beta_ttm(:,m), 'k','LineWidth',1.5)%, i, Beta_LB(:,m), 'b:', i, Beta_UB(:,m),'r:','LineWidth',1.5)
    legend('Time-varying Parameters')%, 'Low Band', 'High Band');
    title('Filtered Time-varying Parameters');


%figure
%plot(i, bm ,'k', i, Beta_ttm, 'b:', 'LineWidth',1.5);
%legend('Time-varying Parameters','Filtered'); 
%title('True and Filtered Time-varying Parameters');

% Smoothed values
    nexttile
    plot(i, Beta_tTm(:,m) ,'k','LineWidth',1.5)%, i, Beta_LB_SM(:,m), 'b:', i, Beta_UB_SM(:,m),'r:','LineWidth',1.5)
    legend('Time-varying Parameters')%, 'Low Band', 'High Band');
    title('Smoothed Time-varying Parameters')% and Confidence Interval');

end
%figure
%plot(i, bm, 'k', i, Beta_tTm,'b:', 'LineWidth',1.5);
%legend('True','Smoothed');
%title('True and Smoothed Time-varying Parameters');


%% State-Space Model
% Dynamic Common Factor Model
% y1(t) = gam1*C(t) + e1(t), e1(t) ~ iidN(0, sig2e1)
% y2(t) = gam2*C(t) + e2(t), e2(t) ~ iidN(0, sig2e2)
% y3(t) = gam3*C(t) + e3(t), e3(t) ~ iidN(0, sig2e3)
% y4(t) = gam4*C(t) + e4(t), e4(t) ~ iidN(0, sig2e4)
% y5(t) = gam5*C(t) + e5(t), e5(t) ~ iidN(0, sig2e5)
% y6(t) = gam6*C(t) + e6(t), e6(t) ~ iidN(0, sig2e6)
% y7(t) = gam7*C(t) + e7(t), e7(t) ~ iidN(0, sig2e7)
% where gam1 = 1

% NO correlation between e1(t), ... ,e7(t)
% C(t) = Mu + phi1*C(t-1) + phi2*C(t-2) + v(t), v(t) ~ iidN(0, sig2v) 

% Measurement equation
% Y(t) = C + H*B(t) + e(t), e(t) ~ iidN(0,R) 

% Transition equation
% B(t) = Mu + F*B(t-1) + v(t), v(t) ~ iidN(0,Q) 

% SS Parameter
% C = [0 0 0 0 0 0 0]'
% H = [1 gam2 gam3 gam4 gam5 gam6 gam7]'
% R = diag([sig2e1; sig2e2; sig2e3; sig2e4; sig2e5; sig2e6; sig2e7])

% Mu = [alpha 0 0]' 
% F = [ phi1 phi2 0 ; 1 0 0 ; 0 1 0 ] 
% Q = diag([sig2v;0;0])

clear; 
clc; 

%% Step 1: DGP
% Data information
T = 1500;     % Sample size 
cut = 500;    % Burn-in size              
k = 1;        % One common factor

% True Parameter
% T.E.
gam2 = 2; 
gam3 = 1; 
gam4 = 2; 
gam5 = 3; 
gam6 = 4; 
gam7 = 1; 

sig2e1 = 0.8; 
sig2e2 = 0.6; 
sig2e3 = 0.4; 
sig2e4 = 0.1; 
sig2e5 = 0.2; 
sig2e6 = 0.3; 
sig2e7 = 0.5; 

% M.E.
mu = 1; 
phi1 = 0.5;
phi2 = 0.3;
sig2v = 0.8;

tru_para = [gam2;gam3;gam4;gam5;gam6;gam7;sig2e1;sig2e2;sig2e3;sig2e4;sig2e5;sig2e6;sig2e7;mu;phi1;phi2;sig2v]; 

% Pre-allocation for c(t)
cm = zeros(T,1);

% Initial c(t), for t = 1
cL = mu + randn(1,1)*sqrt(sig2v);
cm(1) = cL;
cm(2) = mu + phi1*cm(1) + sqrt(sig2v)*randn(1,1);

% Generate c(t) for t = 2 ~ T
for t = 3:T 
    cm(t) = mu + phi1*cm(t-1) + phi2*cm(t-2) + sqrt(sig2v)*randn(1,1);      
end

% Generate y(t)
y1 = cm + sqrt(sig2e1)*randn(t,1); 
y2 2*cm + sqrt(sig2e2)*randn(t,1); 
y3 = gam3*cm + sqrt(sig2e3)*randn(t,1);
y4 = gam4*cm + sqrt(sig2e4)*randn(t,1);
y5 = gam5*cm + sqrt(sig2e5)*randn(t,1);
y6 = gam6*cm + sqrt(sig2e6)*randn(t,1);
y7 = gam7*cm + sqrt(sig2e7)*randn(t,1);= gam

ym = [y1 y2 y3 y4 y5 y6 y7]; 

% Burn-in
ym = ym(cut+1:end,:);
cm = cm(cut+1:end);

%% Step 2: Maxmimum Likelihood Estimation
% Data
data = ym;

% SS Parameter
% C = [0 0 0 0 0 0 0]'
% H = [1 gam2 gam3 gam4 gam5 gam6 gam7]'
% R = diag([sig2e1; sig2e2; sig2e3; sig2e4; sig2e5; sig2e6; sig2e7])

% Mu = [alpha 0 0]' 
% F = [ phi1 phi2 0 ; 1 0 0 ; 0 1 0 ] 
% Q = diag([sig2v;0;0])

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
psi0 = [gam2;gam3;gam4;gam5;gam6;gam7;sig2e1;sig2e2;sig2e3;sig2e4;sig2e5;sig2e6;sig2e7;mu;phi1;phi2;sig2v]; 

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








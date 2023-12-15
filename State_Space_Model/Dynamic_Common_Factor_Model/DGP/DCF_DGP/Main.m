 %% State-Space Model
% Dynamic Common Factor Model
% y(1t) = C(t) + e(1t), e(1t) ~ iidN(0, sig2e1)
% y(2t) = gam*C(t) + e(2t) , e(2t) ~ iidN(0, sig2e2)  
% NO correlation between e(1t) and e(2t)
% C(t) = mu + phi*C(t-1) + v(t), v(t) ~ iidN(0, sig2v) 

% Measurement equation
% y(t) = C + H*B(t) + e(t), e(t) ~ iidN(0,R) 

% Transition equation
% B(t) = Mu + F*B(t-1) + v(t), v(t) ~ iidN(0,Q) 

% SS Parameter
% C = [0 0]', 
% H = [1 gam]',
% R = [sig2e1 0;0 sig2e2]

% Mu = mu, 
% F = phi, 
% Q = sig2v

clear; 
clc; 

%% Step 1: DGP
% Data information
T = 1500;     % Sample size 
cut = 500;    % Burn-in size              
k = 1;        % One common factor

% True Parameter
% T.E.
gam = 3; 
sig2e1 = 2; 
sig2e2 = 1; 

% M.E.
mu = 4; 
phi = 0.5; 
sig2v=0.6;
tru_para = [gam;sig2e1;sig2e2;mu;phi;sig2v]; 

% Pre-allocation for c(t)
cm = zeros(T,1);

% Initial c(t), for t = 1
cL = mu + randn(1,1)*sqrt(sig2v);
cm(1) = cL;

% Generate c(t) for t = 2 ~ T
for t = 2:T 
    cm(t) = mu + phi*cm(t-1) + sqrt(sig2v)*randn(1,1);      
end

% Generate y(t)
y1 = cm + sqrt(sig2e1)*randn(t,1); 
y2 = gam*cm + sqrt(sig2e2)*randn(t,1); 
ym = [y1 y2]; 

% Burn-in
ym = ym(cut+1:end,:);
cm = cm(cut+1:end);

%% Step 2: Maxmimum Likelihood Estimation
% Data
data = ym;

% SS Parameter
% C = 0, 
% H= [1 gam]',
% R = [sig2e1 0;0 sig2e2]

% Mu = mu, 
% F = phi, 
% Q = sig2v

% Block for each parameters
indH = 1;
indR = [2;3];
indMu = 4;
indF = 5;
indQ = 6;

% Structure variables
Sn.data = data;
Sn.indH = indH;
Sn.indR = indR;
Sn.indMu = indMu;
Sn.indF = indF;
Sn.indQ = indQ;

% Initial values
psi0 = [gam;log(sig2e1);log(sig2e2);mu;phi;log(sig2v)]; 

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








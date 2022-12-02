%% Maximum likelihood estimation
% Autoregressive Model(2)
% Y(t) = mu + phi*Y(t-1) + e(t), e(t)~iidN(0,sig2)

clear;
clc;
addpath('C:\Users\EWHA\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: DGP
% data load
Ym = xlsread('Kospi.csv','Kospi','C3:C5560');
Zm = xlsread('KOSPI.xlsx', '데이터', 'B2:B132');

%% Step 2: Estimation 
Y = Ym(4:end);    % T by 1
YL1 = Ym(3:end-1);
YL2 = Ym(2:end-2);
YL3 = Ym(1:end-3);

T = rows(Y);
X = [ones(T,1) YL1 YL2 YL3];      % T by 2
k = cols(X);

% Estimation
% parameter (c, phi, sig2)
theta0 = [0;0.1;0.1;0.1;0.1];
Data = [Y X];

% index
index = [1;2;3;4;5];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);

parcorr(Ym);
rhoja = autocorr(Y, T-1);
rhojb = autocorrr(Y, T-2);

% Compute t-value and p-value
diag_cov = diag(Vj);
stde = sqrt(diag_cov);
t_val = thetamx./stde; 
p_val = 2*(1 - cdf('t', abs(t_val), T-k)); 

% Compute [AIC BSC R2]
mY = meanc(Y); % sample mean of the dependent variable
bhat = thetamx(1:4,1);
Yhat = X*bhat; % T by 1, fitted value
ehat = Y - Yhat; % T by 1, residuals
TSS = Y'*Y - T*mY^2; % TSS
RSS = ehat'*ehat;             % RSS
AIC = log(RSS/T) - 2*k/T;     % AIC
BSC = log(RSS/T) + k*log(T)/T; % BSC
R2_ = 1 - (T-1)*RSS/(TSS*(T-k)); % Adjusted R2


[h1, pValueadf] = adftest(Y);
[h2, pvaluekpss] = kpsstest(Y);
[h3, pValuelbq] = lbqtest(ehat, 'lags', T-1);

[h4, pValueadf1] = adftest(Zm);
[h5, pValuekpss1] = adftest(Zm);
%% Step 3: Results %%
disp('=========================================');
disp(['  Estimates ','  t value ', '  p value ']); 
disp('=========================================');
disp([thetamx t_val p_val]); 
disp('=========================================');
disp('     ADF       KPSS    Ljung-Box');
disp('=========================================');
disp([pValueadf, pvaluekpss, pValuelbq]);
disp('=========================================');
disp('    AIC       BSC          R2_');
disp('=========================================');
disp([AIC BSC R2_]);
disp('=========================================');
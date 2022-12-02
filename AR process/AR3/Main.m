%% Maximum likelihood estimation
% Autoregressive Model(2)
% Y(t) = mu + phi*Y(t-1) + e(t), e(t)~iidN(0,sig2)

clear;
clc;
addpath('C:\Users\EWHA\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: DGP
% data load
Ym = xlsread('KOSPI','데이터','B2:B133');

%% Step 2: Estimation 
Y = Ym(4:end);             % T by 1
YL1 = Ym(3:end-1);
YL2 = Ym(2:end-2);
YL3 = Ym(1:end-3);
T = rows(Y);
X = [ones(T,1) YL1, YL2, YL3];      % T by 2
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

% Compute t-value and p-value
diag_cov = diag(Vj);
stde = sqrt(diag_cov);
t_val = thetamx./stde; 
p_val = 2*(1 - cdf('t', abs(t_val), T-k)); 

%% Step 3: Results %%
disp('=========================================');
disp(['  Estimates ','  t value ', '  p value ']); 
disp('=========================================');
disp([thetamx t_val p_val]); 
disp('=========================================');

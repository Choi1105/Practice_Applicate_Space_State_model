%% Maximum likelihood estimation
% Autoregressive Model(2)
% Y(t) = mu + phi1*Y(t-1)  phi2*Y(t-2)+ e(t), e(t)~iidN(0,sig2)

clear;
clc;
addpath('C:\Users\KYM\Documents\MATLAB\AR process\AR2');

%% Step 1: DGP
% data load
Ym = readmatrix("spinach0.csv",Range="B2:B7668");

%% Step 2: Estimation 
Y = Ym(3:end);             % T by 1
YL1 = Ym(2:end-1);
YL2 = Ym(1:end-2);
T = rows(Y);
X = [ones(T,1) YL1, YL2];      % T by 2
k = cols(X);

% Estimation
% parameter (c, phi1,phi2 sig2)
theta0 = [0;0.1;0.1;0.1];
Data = [Y X];

% index
index = [1;2;3;4];
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

%% Maximum likelihood estimation
% Autoregressive Model(1)
% Y(t) = mu + phi*Y(t-1) + e(t), e(t)~iidN(0,sig2)

clear;
clc;
addpath('C:\Users\KYM\Documents\MATLAB\AR process\AR1');

%% Step 1: DGP
% data load
Ym = readmatrix("spinach0.csv",Range="B2:B7668");

%% Step 2: Estimation 
Y = Ym(2:end);             % T by 1
YL = Ym(1:end-1);
T = rows(Y);
X = [ones(T,1) YL];      % T by 2
k = cols(X);
rhom = zeros(T,1);

% Estimation
% parameter (c, phi, sig2)
theta0 = [0;0.5;0.1];
Data = [Y X];

% index
index = [1;2;3];
printi = 1;

% Optimization

[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);

rhoj = autocorr(Y, 128);
[y, lags] = xcorr(Y, Y);

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

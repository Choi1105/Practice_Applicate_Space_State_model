%% Dynamic Forecasting
% Iterative
clear;
clc;

%% Step 1: Data Load
Data = xlsread('project','sheet1','B2:B260');

% Data
Y_full = Data(2:end,1);

X1 = ones(rows(Y_full),1);
X2 = Data(1:end-1,1);
X_full = [X1 X2];

%% Step 2: Estimation
Y = Y_full(1:end-1);
X = X_full(1:end-1,:);
printi = 0;
[bhat, Yhat, ehat, sig2hat, varbhat, stde, t_val, TSS, RSS, R2, R2_, SC, AIC] = OLS_OUT(Y,X,printi)

%% Step 3: Dynamic Forecasting
% Actual
Ya = Y_full(end);

% Forecasting
Xf = X_full(end,:)';
Yf_11_Nov = Xf'*bhat;
Yf_12_Nov = [1 Yf_11_Nov]*bhat;
Yf_15_Nov = [1 Yf_12_Nov]*bhat;
Yf_16_Nov = [1 Yf_15_Nov]*bhat;
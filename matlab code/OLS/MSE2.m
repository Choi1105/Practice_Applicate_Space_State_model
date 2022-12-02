%% Static Forecasting

clear;
clc;

%% Step 1: Data Load, 2016-11-14 - 2021-11-11
Data = xlsread('project','sheet1','B2:B260');

% Data
Y_full = Data(2:end,1);

X1 = ones(rows(Y_full),1);
X2 = Data(1:end-1,1);
X_full = [X1 X2];

%% Step 2: 2021-11-10 Forecasting
Y = Y_full(1:end-2);
X = X_full(1:end-2,:);
printi = 0;
[bhat, Yhat, ehat, sig2hat, varbhat, stde, t_val, TSS, RSS, R2, R2_, SC, AIC] = OLS_OUT(Y,X,printi)

% Actual
Ya = Y_full(end-1);

% Forecasting
Xf = X_full(end-1,:)';
Yf = Xf'*bhat;

% Forecast Error
FE1 = (Ya - Yf);

%% Step 3: 2021-11-11 Forecasting
Y = Y_full(1:end-1);
X = X_full(1:end-1,:);
printi = 0;
[bhat, Yhat, ehat, sig2hat, varbhat, stde, t_val, TSS, RSS, R2, R2_, SC, AIC] = OLS_OUT(Y,X,printi)

% Actual
Ya = Y_full(end);

% Forecasting
Xf = X_full(end,:)';
Yf = Xf'*bhat;
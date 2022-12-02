%% Linear Regression Model
% OLS

clear
clc

%% Step 1: Data Generating Process
T = 500;
b1 = 1;
b2 = 2;
b3 = 3;
sig2 = 1;

X1m = ones(T,1);
X2m = 5*rand(T,1);
X3m = 5*rand(T,1);
em = sqrt(sig2)*randn(T,1);

Ym = b1*X1m + b2*X2m + b3*X3m + em;

%% Step 2: Estimation
% Data
Y = Ym;
X = [X1m X2m X3m];

% OLS estimator
bhat = inv(X'*X)*X'*Y;
Yhat = X*bhat
ehat = Y - Yhat;

% Figure
plot([Y Yhat]);
legend('Y','Yhat');
title('Y and Yhat');
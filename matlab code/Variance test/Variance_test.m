%% Variance test

clear
clc

%% Step 1: Data Generating
T = 40;
B0 = 1;
B1 = -1;
sig2 = 1;

X1m = ones(T,1);
X2m = randn(T,1);
em = sqrt(sig2)*randn(T,1);
B11m = zeros(5000,1);
B12m = zeros(5000,1);

Ym = X1m*B0 + X2m*B1 + em;

% Data
Y = Ym;
X = [X1m X2m];

% OLS estimator
bhat1 = inv(X'*X)*X'*Y;
bhat2 = Y(2,1)-Y(1,1)/X2m(2,1)-X2m(1,1);
Yhat = X*bhat1;
ehat = Y - Yhat;

% Figure
plot([Y Yhat]);
legend('Y','Yhat');
title('Y and Yhat');

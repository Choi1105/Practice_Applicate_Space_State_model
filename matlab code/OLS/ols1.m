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

%% Step 2 : Estimation
% Data
Y = Ym;
X = [X1m X2m X3m];
T = rows(X);
k = cols(X);

% OLS estimator
bhat = inv(X'*X)*X'*Y;
Yhat = X*bhat;
ehat = Y - Yhat;
sig2hat = ehat'*ehat/(T-k);
varbhat = sig2hat*inv(X'*X);
stde = sqrt(diag(varbhat));

% Model comparison
RSS = ehat'*ehat;
TSS = (Y - mean(Y))'*(Y - mean(Y));
R2 = 1 - RSS/TSS
R2_ = 1 - (RSS*(T-1))/(TSS*(T-k));
SC = log(RSS/T) + k*log(T)/T;
AIC = log(RSS/T) + 2*k/T;

% Table
disp('========================');
disp([' Estimates ', '  S.E  '])
disp('========================');
disp([bhat stde]);
disp('========================');
disp(['S.E. of regression is ', num2str(sqrt(sig2hat))]);
disp(['R2 is ', num2str(R2)]);
disp(['adjusted R2 is ', num2str(R2_)]);
disp(['SC is ', num2str(SC)]);
disp(['AIC is ', num2str(AIC)]);
disp('========================');

% Figure
plot ([Y Yhat]);
legend('Y', 'Yhat');
title('Y and Yhat');
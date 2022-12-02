%% Gibbs Sampling
% Multi-Step ahead Forecasting based on Linear Regression model

%% Step 0: Set up
clear;
clc;
addpath('C:\Users\이건희\OneDrive - 전북대학교\문서\MATLAB\M_library\M_library');

%% Step 0 : Data Generating Process
% Sample size
T = 200;

% True value
beta1 = 1;
beta2 = 2;
sig2 = 0.5;

% Independent/Dependent variable
X1m = ones(T,1);
X2m = 5*rand(T,1);

Ym = X1m*beta1 + X2m*beta2 + sqrt(sig2)*randn(T,1);

%% Step 1: Prior
% Full Sample
Y_Full = Ym;                % T by 1
X_Full = [X1m, X2m];        % T by 2  

% Information
[T, k] = size(X_Full);

% beta prior
beta_0 = zeros(k,1);
B_0 = 100*eye(k);

% sig2 prior
% mean = d_0/v_0
v_0 = 3;
d_0 = 3;

%% Step 2: Sampling
n0 = 5000;
n1 = 5000;
H = 25 ;
[betam,sig2m,yfm,SEm,PPD] = Bayes_linear_N_Pred(Y_Full,X_Full,beta_0,B_0,v_0,d_0,n0,n1,H);

%% Step 3 : Results
RMSE = sqrt(meanc(SEm));
lnPPL = sumc(log(PPD));
disp(['RMSE = ', num2str(RMSE)]);
disp(['Log PPL = ',num2str(lnPPL)]);

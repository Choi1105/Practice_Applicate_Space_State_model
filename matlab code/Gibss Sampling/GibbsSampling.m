%% Gibbs Sampling
% Linear Regression model
% Y = b1 + x2*b2 + x3*b3 + e, Y = Xb + e, e ~ iidN(0, Sig2)

%% Step 0: Set up
clear;
clc;
addpath('C:\Program Files\MATLAB\R2021b\M_library');

% Generate Data
T = 200;

% True value
beta1 = 1;
beta2 = 2;
beta3 = 3;
sig2 = 0.5;

X1m = ones(T,1);
X2m = 5*rand(T,1);
X3m = 3*rand(T,1);

Ym = X1m*beta1 + X2m*beta2 + X3m*beta3 +sqrt(sig2)*randn(T,1);

%% Step 1: Prior
% Data
Y = Ym;                     % T by 1
X = [X1m, X2m, X3m];        % T by 3  

% Information
[T, k] = size(X);

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
[betam,sig2m] = Bayes_linear_N(Y,X,beta_0,B_0,v_0,d_0,n0,n1);

%% Step 3: Result
MHm = [betam sig2m];

alpha = 0.05;
maxac = 200;
is_plot = 0;
postmom = MHout(MHm, alpha, maxac, is_plot);

% Display
disp('==========================================');
disp(['   Mean''   S.E.' '     2.5%' '     97.5%']);
disp('==========================================');
disp([postmom(:,2) postmom(:,3) postmom(:,4) postmom(:,6)]);


%% Step 4: OLS Estimation and Results
printi = 1;
[bhat, sig2hat, stde, t_val, Yhat, ehat, varbhat, mY, TSS, RSS, R2, R2_,...
    SC, AIC] = OLSout(Y,X,printi);
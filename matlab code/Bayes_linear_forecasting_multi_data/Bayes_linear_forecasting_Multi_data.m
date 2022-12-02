%% Gibbs Sampling
% Multi-Step ahead Forecasting based on Linear Regression model

%% Step 0: Set up
clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 0 : Data Generating Process
ym = xlsread('Oil price data', 'Oil price', 'K230:K408');  % Oil Price (t)

T = rows(ym);

x1m = ones(T,1);
x2m = xlsread('원유생산및소비', 'sheet1', 'C2:C180');        % Oil Production(t)
x3m = xlsread('Oil price data', 'Oil price', 'K229:K407');  % Oil Price (t)
x4m = xlsread('Oil price data', 'Oil price', 'M229:M407');  % Oil Consumption (World Industrial Price Index) (t)

%% Step 1: Prior
% Full Sample
Y = ym;
X = [x1m, x2m, x3m, x4m];

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
H = 1;
[betam,sig2m,yfm,SEm,PPD] = Bayes_linear_N_Pred(Y,X,beta_0,B_0,v_0,d_0,n0,n1,H);

%% Step 3 : Results
RMSE = sqrt(meanc(SEm));
lnPPL = sumc(log(PPD));
disp(['RMSE = ', num2str(RMSE)]);
disp(['Log PPL = ',num2str(lnPPL)]);

%% Gibbs Sampling
% Linear Regression model
% Y = b1 + x2*b2 + e, Y = Xb + e, e ~ iidN(0, Sig2)

%% Step 0: Set up
clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library\M_library');

% Data Load
Ym = readmatrix('KOSPI','sheet','데이터','range','B2:B110');              % KOSPI 지수 변화율
X1m = ones(109,1);
X2m = readmatrix('KOSPI','sheet','데이터','range','C2:C110');             % PMI 변화율
X3m = readmatrix('KOSPI','sheet','데이터','range','D2:D110');             % 환율 변화율
 
%% Step 1: Prior
% Full Sample
Y_Full = Ym;                     % T by 1
X_Full = [X1m, X2m, X3m];        % T by 3  

% In-sample data                 % 예측하기 전의 data
Y = Y_Full(1:end-1);
X = X_Full(1:end-1,:);

% Information
[T, k] = size(X);

% beta prior
beta_0 = zeros(k,1);
B_0 = 1*eye(k);

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

%% Step 4: Forecasting
ya = Y_Full(end);     % Actual
xf = X_Full(end,:)';  % Exogenous variables are given

% Forecasting
yfm = zeros(n1,1);
for iter = 1:n1
    beta_i = betam(iter,:)';
    sig2_i = sig2m(iter);
    yfm(iter) = mvnrnd(xf'*beta_i, sig2_i);
end

%% Step 5: Evaluation
% SE
SE = (ya - meanc(yfm))^2;

% PPD
PPDm = zeros(n1,1);
for iter = 1:n1
    beta_i = betam(iter,:)';
    sig2_i = sig2m(iter);
    PPDm(iter) = normpdf(ya, xf'*beta_i, sqrt(sig2_i));
end
PPD = meanc(PPDm);

% Display
disp(['SE = ',num2str(SE)]);
disp(['PPD = ',num2str(PPD)]);


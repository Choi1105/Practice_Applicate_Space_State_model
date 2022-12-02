%% Recursive BVAR model
% Short-Run restrictions based on Cholesky Decomposition
% Simulation, SVAR(2)

clear;
clc;
addpath('C:\Users\이건희\OneDrive - 전북대학교\문서\MATLAB\M_library')

%% Step 1: Data
% Data information
y1m = xlsread('변동성측정', 'Sheet1', 'E2:E421');  % Oil Production (t)  
y2m = xlsread('변동성측정', 'Sheet1', 'F2:F421');  % Oil Price (t)
y3m = xlsread('변동성측정', 'Sheet1', 'G2:G421');  % Oil Consumption (World Industrial Price Index) (t)

%% Step 1: Prior
% Full Sample
Data = [y1m y2m y3m];

%% Step 2: Set-up
% Set VAR order: p
% Max lag for variance decomposition (VD) and impulse response (IR)
p = 1;
mlag = 60;

% Number of variables
k = cols(Data);

% Demean Data
Y = demeanc(Data);

% Burn in(n0) / Actual sampling size (n1)
n0 = 1000;
n1 = 2000;
n = n0 + n1;

% Conjugate Prior Distributions
% Normal prior for beta(-vec(phi))
pkk = p*k*k;
beta0 = zeros(pkk,1);
B0 = 0.1*eye(pkk);

% Inverse Wishart prior for inv(Omega)
% Mean = R0/(nu - p - 1)
nu0 = 20;
R0 = 0.2*eye(k);

%% Step 3: Gibbs Sampling
[PSIm, MHm] = Recursive_VAR(n0, n1, beta0, B0, nu0, R0, Y, p, mlag);

%% Step 4: Report Posterior Distributions
% criteria for output
alpha = 0.025;        % 신뢰구간/2
maxac = 200;          % autocorrelation의 최대 숫자
is_postplot = 0;      % -1, display posterior distribution, -0, NOT

% Impulse response function
Plot_IRF(PSIm, alpha);

% Posterior Moments / Sampling efficiency
postmom = MHout(MHm, alpha, maxac, is_postplot);

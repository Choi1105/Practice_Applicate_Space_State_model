%% Recursive BVAR model
% Short-run restrictions based on Cholesky Decomposition
% Simulation, SVAR(2)

clear;
clc;
addpath('C:\Users\이건희\OneDrive - 전북대학교\문서\MATLAB\M_library')

%% STEP 1: Data Generating Process
% Data information
% Actual sample size: T - cut
T = 800;
cut = 500;
k = 2;
p = 2;

% True parameters
phi1 = 0.8*eye(k);
phi2 = 0.15*eye(k);
Omega = 0.5*eye(k);
% Data Generating, VAR(p)
ym = zeros(T,k);
for t = p+1:T
yt = phi1*ym(t-1,:)' + phi2*ym(t-2,:)' + chol(Omega)'*randn(k,1);
ym(t,:) = yt';
end
% Burn-in initial data
Data = ym(cut:end,:);

%% Step 2: Set-up
% Set VAR order: p
% Max lag for variance decomposition (VD) and impulse response (IR)
p = 2;
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

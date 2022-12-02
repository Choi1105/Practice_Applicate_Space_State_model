clear
clc
addpath('C:\Users\이건희\OneDrive - 전북대학교\문서\MATLAB\M_library')

%% Step 1: Data Generating Process
% information
T = 300;
k = 200;

% parameters
beta = zeros(k,1);
onesk = 10:10:k;
beta(onesk) = 10;
sig2 = 1;

% data
Xm = 5*randn(T, k);
Ym = Xm*beta + sqrt(sig2)*randn(T,1);

%% Step 1: Prior
Y = Ym;
X = Xm;
k = cols(X);

% sig2 prior
v_0 = 3;
d_0 = 3;

% p prior
c0 = 5;
d0 = 5;

% b0, b1 prior  가장 중요한 prior -> b0와 b1의 크기에 따라 중요한 변수와 중요하지 않은 변수를 선택하는 기준
v_00 = 10;             % prior의 parameter를 변화시키는 것도 bayesian 분석에서 새로운 모형을 적용시키는 것
d_00 = (0.01^2)*10;    % 최적의 prior를 찾는 것도 중요한 연구 중 하나

v_01 = 10;
d_01 = 1*10;

%% Step 2: Sampling
n0 = 1000;    % burn-in
n1 = 1000;    % MCMC size
[betam, sig2m, Gam, postmom_beta, postmom_sig2, postmom_Gam] = MCMC_SSVS(Y, X, n0, n1, v_0, d_0, v_00, d_00, v_01, d_01, c0, d0);

%% Step 3: Results
Gam_hat = postmom_Gam(:, 2);
figure
x = 1:rows(Gam_hat);
scatter(x, Gam_hat, 50, 'b', 'filled')
ylim([-0.05, 1.05])
xlim([0, length(x)+1])
xlabel('Variable')
ylabel('Prob. of inclusion')
grid on
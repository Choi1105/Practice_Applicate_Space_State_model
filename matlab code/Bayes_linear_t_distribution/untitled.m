%% Gibbs Sampling
% 오차항이 스튜던트 t분포인 모형
% beta ~ N(beta0, B0)
% sig2 ~ IG(a_0/2, d_0/2)
% lamda(t) ~ Gamma(v/2, v/2)
% Y = b1 + x2*b2 + e, Y = Xb + e, e ~ N(0, Sig2)

%% Step 0: Set up
clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: DGP
T = 100;

% True value
beta1 = 1;
beta2 = 2;
sig2 = 0.5;

% Independent/Dependent variable
x1m = ones(T,1);
x2m = 5*rand(T,1);
nu = 10;

ym = x1m*beta1 + x2m*beta2 + sqrt(sig2)*randn(T,1);

%% Step 2: Prior
Y = ym;
X = [x1m, x2m];
Beta = [beta1, beta2];

Y = zeros(T,1);
for t = 1:T
    xt = X(t,:)';
    lam = randgam(nu/2, nu/2);
    Y(t) = xt'*Beta + sqrt(sig2/lam)*randn(1,1);
end
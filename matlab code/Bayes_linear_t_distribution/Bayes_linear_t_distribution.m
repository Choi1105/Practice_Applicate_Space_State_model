%% Gibbs Sampling
% 오차항이 스튜던트 t분포인 모형
% beta ~ N(beta0, B0)
% sig2 ~ IG(a_0/2, d_0/2)
% lamda(t) ~ Gamma(v/2, v/2)
% Y = X*Beta + e, Y = X*Beta + e, e ~ N(0, Sig2)

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
X = [x1m, x2m];
Beta = [beta1; beta2];

Y = zeros(T,1);
for t = 1:T
    xt = X(t,:)';
    lambda = randgam(nu/2, nu/2);
    Y(t) = xt'*Beta + sqrt(sig2/lambda)*randn(1,1);
end

%% Step 2: Prior
k = cols(X);
beta_0 = zeros(k,1);
B_0 = 100*eye(k);

alpha_0 = 10;
delta_0 = 10;

n0 = 1000;
n1 = 1000;
n = n0 + n1;

sig2 = delta_0/alpha_0;
Lambda = eye(T);

betam = zeros(n,k);
sig2m = zeros(n,1);

for iter = 1:n 

    %% beta Sampling
    B1 = inv((1/sig2)*X'*lambda*X + inv(B_0));
    A = (1/sig2)*X'*lambda*Y + inv(B_0)*beta_0;
    beta = mvnrnd(B1*A, B1)';
    betam(iter,:) = beta';
   
    %% sig2 Sampling
    alpha_1 = alpha_0 + T;
    delta_1 = delta_0 +(Y - X*Beta)'*Lambda*(Y - X*Beta);
    sig2 = randig(alpha_1/2, delta_1/2);
    sig2m(iter,:) = sig2';

    %% lamda Sampling
    for t = 1:T
        lambda = randgam((nu+1)/2, 1/2*(nu + (1/sig2)*(Y(t) - X(t,:)*Beta)^2));
        Lambda(t,t) = lambda;

    end

end

% Burn-in
betam = betam(n0+1:n,:);
sig2m = sig2m(n0+1:n,1);

MHm = [betam sig2m];
alpha = 0.05;
maxac = 200;
is_plot = 1;
postmom = MHout(MHm, alpha, maxac, is_plot);

% Display
disp('==========================================');
disp(['   Mean''   S.E.' '     2.5%' '     97.5%']);
disp('==========================================');
disp([postmom(:,2) postmom(:,3) postmom(:,4) postmom(:,6)]);

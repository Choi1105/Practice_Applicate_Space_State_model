%% Markov-switching Model
% 2 State : at = 1, 2
% y(t) = b1_s(t) + x2(t)*b2_s(t) + e(t)
% where, e(t) ~ i.i.d.N(0, sig2_s(t))


clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: Data load
ym = readmatrix('한전데이터','sheet','Sheet1','range','I4:I67');
x1m = ones(64,1);
x2m = readmatrix('한전데이터','sheet','Sheet1','range','J4:J67');
x3m = readmatrix('한전데이터','sheet','Sheet1','range','K4:K67'); 
x4m = readmatrix('한전데이터','sheet','Sheet1','range','L4:L67');
x5m = readmatrix('한전데이터','sheet','Sheet1','range','M4:M67'); 

%% Step 1: Prior
% Full Sample
Y = ym;       % 영업이익 차분    
X = [x1m, x2m, x3m, x4m, x5m];

T = rows(Y);  % Sample size
k = cols(X);  % 변수 개수
ns = 2;       % regime 개수

sig21 = 0.1;
sig22 = 0.5;


% S(0)
sm = zeros(T,1);

for t = 1:T
if rand(1,1) < 0.5; st = 1;

else                st = 2;

end

sm(t,1) = st;
end


% Sampling size (n0: burn-in, n1: actual sampling size)
n0 = 1000;
n1 = 1000;
n = n0 + n1;

% Prior for beta, Normal
beta0 = zeros(k, 1);
B0 = 1*eye(k);

% Prior for beta, Inverse Gamma
v_0 = 2;
d_0 = 2;

% Prior for P, Beta
% p11
a10 = 10; b10 = 10;
% p22
a20 = 10; b20 = 10;

% Initial values
S = sm;
p11 = a10/(a10 + b10);
p22 = a20/(a20 + b20);
sig21L = sig21;
sig22L = sig22;

P = zeros(ns,ns);
P(1,1) = p11;
P(1,2) = 1 - p11;
P(2,2) = p22;
P(2,1) = 1 - p22;

% Pre-allocation for parameters
Sm = zeros(n, T);
beta1m = zeros(n, k);
beta2m = zeros(n, k);
sig2m = zeros(n, ns);
Pm = zeros(n, ns);

%% Step 3: MCMC Sampling
for iter = 1:n

    % Stpe 3-1: Separate samples
    % When s(t) = 1
    Y1 = Y(S == 1);
    X1 = X(S == 1, :);

    % When s(t) = 2
    Y2 = Y(S == 2);
    X2 = X(S == 2, :);

    % Step 3-2: beta sampling
    beta1 = Gen_beta(Y1, X1, beta0, B0, sig21);
    beta2 = Gen_beta(Y2, X2, beta0, B0, sig22);
    beta1m(iter, :) = beta1';
    beta2m(iter, :) = beta2';

    % Step 3-3: sig2 sampling
    sig21 = Gen_sig2(Y1, X1, v_0, d_0, beta1);
    sig22 = Gen_sig2(Y2, X2, v_0, d_0, beta2);

    % Identification restrictions
    if sig21 > sig22
        sig21 = sig21L;
        sig22 = sig22L;
    end
    sig2m(iter, :) = [sig21, sig22];

    % Updating
    sig21L = sig21;
    sig22L = sig22;

    % Step 3-4: S sampling
    Filtpm = Hamilton_Filter(Y, X, beta1, beta2, sig21, sig22, p11, p22);
    S = sgen(Filtpm, P);
    Sm(iter, :) = S';

    % Step 3-5: P sampling
    [p11, p22, P] = Gen_Prb(S, a10, b10, a20, b20);
    Pm(iter, :) = [p11, p22];

end

%% Stpe 4: Results
% Burn-in
MHm = [beta1m, beta2m, sig2m, Pm];
MHm = MHm(n0+1:n,:);

% Summary posterior distributions
is_postplot = 1;
alpha = 0.025;
maxac = 200;
postmom = MHout(MHm, alpha, maxac, 0);

% Display Figure
Sm = Sm(n0+1:n, :);
post_prob2 = meanc(Sm) - 1;
i = 1:T;
i = i';
hold on;
plot(i, (sm-1),  'LineStyle', '-', 'LineWidth', 2.0, 'Color', 'b');
plot(i, post_prob2,  'LineStyle', '--', 'LineWidth', 1.5, 'Color', 'r');
plot(i, 1.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
plot(i, -0.01*ones(T,1),  'LineStyle', ':', 'LineWidth', 1.0, 'Color', 'k');
legend('DGP', 'Estimates')
title('Posterior Prob. of Regime 2');
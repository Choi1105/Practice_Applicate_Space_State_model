%% Probit model
% beta ~ N(beta0, B0)
% z(t) ~ N(xt*beta, 1)
% y(t) = I (z(t) > 0)


%% Step 0: Set up
clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');


%% Step 1: DGP
T = 500;
x1 = ones(T,1);
x2 = 5*randn(T,1);
x3 = 5*randn(T,1);

X = [x1, x2, x3];

% True parameter
beta1 = -0.5;
beta2 = -2;
beta3 = 3;

beta = [beta1; beta2; beta3];

Z = X*beta + randn(T, 1);

Y = Z > 0; % Z가 양수면 1, 음수면 0

%% Step 2: Prior
k = cols(X);

beta0 = zeros(k,1);
B0 = eye(k);
B0inv = inv(B0);

n0 = 200000;
n1 = 200000;
n = n0 + n1;

Z = demeanc(Y);  % Z의 초기값, demeanc = 평균을 제거
betam = zeros(n, k);

for iter = 1:n
    %% Beta Sampling
    B1 = inv(X'*X + B0inv);  % sig2 = 1
    A = X'*Z + B0inv*beta0;
    beta = mvnrnd(B1*A, B1)';
    betam(iter,:) = beta'; 

    %% Z Sampling
    Z = Y.*trandn(X*beta, 1, 0, 20) + (1 - Y).*trandn(X*beta, 1, -20, 0); % Y가 1일때는 Z가 양수, 0일때는 Z가 음수!

end

% Burn-in
betam = betam(n0+1:n,:);


MHm = [betam];
alpha = 0.05;
maxac = 200;
is_plot = 1;
postmom = MHout(MHm, alpha, maxac, is_plot);

% Display
disp('==========================================');
disp(['   Mean''   S.E.' '     2.5%' '     97.5%']);
disp('==========================================');
disp([postmom(:,2) postmom(:,3) postmom(:,4) postmom(:,6)]);

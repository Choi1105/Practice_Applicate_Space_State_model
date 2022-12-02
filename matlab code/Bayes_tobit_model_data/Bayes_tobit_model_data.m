%% tobit model
% beta ~ N(beta0, B0)
% z(t) ~ N(xt*beta, 1)
% y(t) = z(t) (y(t) > 0)
% y(t) = 0 (z(t) < 0)


%% Step 0: Set up
clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');


%% DGP
Data = xlsread('0519.xlsx','2019','G2:AW2289');

income = Data(:,1);
assets = Data(:,2);
wage = Data(:,4);
age = Data(:,7);
d_digital = Data(:,15);
benke = Data(:,11);
%% yuce
Y = benke;
T = rows(Y);
X = [ones(T,1) income wage age];
Z = demeanc(Y);  % Z의 초기값, demeanc = 평균을 제거

Y = zeros(T, 1);
for t = 1:T
    if Z(t) > 0
        Y(t) = Z(t);
    end

end

%% Step 2: Prior
[T,k] = size(X);

beta0 = zeros(k,1);
B0 = eye(k);
B0inv = inv(B0);

a_0 = 3;
d_0 = 3;

sig2 = a_0/d_0;

n0 = 10000;
n1 = 10000;
n = n0 + n1;

betam = zeros(n, k);
sig2m = zeros(n, 1);

%% Step 3: Sampling

for iter = 1:n
    %% Beta Sampling
    B1 = inv(1/sig2*X'*X + B0inv);  % sig2 = 1
    A = 1/sig2*X'*Z + B0inv*beta0;
    Beta = mvnrnd(B1*A, B1)';
    betam(iter,:) = Beta'; 

    %% sig2 Sampling
    a_1 = a_0 + rows(Y);
    d_1 = d_0 + (Z-X*Beta)'*(Z-X*Beta);
    sig2 = randig(a_1/2,d_1/2,1,1);
    sig2m(iter,1) = sig2;

    %% Z Sampling
    for t = 1:T
        xt = X(t, :)';
        mu = xt'*Beta;

        if Y(t) > 0
            Z(t) = Y(t);
        
        else Y(t) == 0;
            Z(t) = trandn(mu, 1, -100, 0);

        end

    end

end

%% Step 4: conclusion?

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

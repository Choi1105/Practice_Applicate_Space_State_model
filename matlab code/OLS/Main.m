%% Linear Regression Model
% OLS
clear
clc

%% Step 1: Data Genereating Process
T = 500;

b1 = 1;
b2 = 2;
b3 = 3;
sig2 = 1;

X1m = ones(T,1);
X2m = 5*rand(T,1);
X3m = 5*rand(T,1);
em = sqrt(sig2)*randn(T,1);

Ym = b1*X1m + b2*X2m + b3*X3m + em;

%% Step 2: Estimation
% Data
Y = Ym;
X = [X1m X2m X3m];
printi = 0;
[bhat, Yhat, ehat, sig2hat, varbhat, stde, t_val, TSS, RSS, R2, R2_, SC, AIC] = OLS_OUT(Y,X,printi);

%% Step 3: Joint Hypothesis
% Test #1 : b1 = b2 = b3 = 0
R = [1 0 0; 0 1 0; 0 0 1];
gam = [0;0;0];

T = rows(X);
k = cols(X);
L = rows(R);

F_val = (R*bhat - gam)'*inv(R*varbhat*R')*(R*bhat - gam)/L;
P_val = 1 - cdf('f', F_val, L, T-k);

% Test #2: b1 + b2 = 1 and b2 = b3
R = [1 1 0; 0 1 -1];
gam = [1; 0];

T = rows(X);
k = cols(X);
L = rows(R);

F_val = (R*bhat - gam)'*inv(R*varbhat*R')*(R*bhat - gam)/L;
p_val = 1 - cdf('f', F_val, L, T-k);

%% Step 4: Results %%
disp('=========================================================')
disp('  Joint  Hypothesis');
disp('=========================================================')
disp(['  H0        ', '                 F_Value', '      P_Value']);
disp('=========================================================')
disp([' b1 = b2 = b3 = 0            ', num2str(F_val), '      ', num2str(P_val)]);
disp('=========================================================')
disp([' b1 + b2 = 1        ', '        b2 = b3     ', num2str(F_val),'        ',num2str(P_val)])

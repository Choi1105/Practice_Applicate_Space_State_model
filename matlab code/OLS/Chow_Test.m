%% Linear Regression Model
% OLS
clear
clc

%% Step 1: Data Loading
Data = xlsread('HW1','VAR','B2:C318');

%% Step 2: Estimation
% Data
Y = Data (:,1);
T = rows(Data);

% Figure
datat = 1:T;
xtick = 13:48:T;
xticklabel = {'Jan 96', 'Jan 00', 'Jan 04', 'Jan 08', 'Jan 12', 'Jan 16', 'Jan 20'};
plot(datat, Y(:,1), 'b-', 'Linewidth',1.5);
set(gca, 'XTick', xtick);
set(gca, 'XTickLabel', xticklabel, 'Fontsize',10);
xlim([0 T+1])

% 구조변화 시점: 2008년 12월
tau = 165;
D1 = [ones(tau,1);zeros(T-tau,1)];
D2 = [zeros(tau,1);ones(T-tau,1)];
X = [D1,D2];
k = cols(X);

printi = 0;
[bhat, Yhat, ehat, sig2hat, varbhat, stde, t_val, TSS, RSS, R2, R2_, SC, AIC] = OLS_OUT(Y,X,printi)

%% Step 3: Joint Hypothesis
% Test #1: b1 = b2
R = [1 -1];
gam = 0;
L = rows(R);

F_val = (R*bhat - gam)'*inv(R*varbhat*R')*(R*bhat - gam)/L;
P_val = 1 - cdf('f', F_val, L, T-k);

%% Step 4: Results
disp('========================================');
disp('Joint Hypothesis');
disp('----------------------------------------');
disp(['  H0     ','           F_Value','       P_Value']);
disp('----------------------------------------');
disp([' Mu1 = Mu2       ', num2str(F_val),   num2str(P_val)]);
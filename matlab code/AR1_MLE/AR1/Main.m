%% Maximum likelihood estimation
% Autoregressive Model(1)
% Y(t) = mu + phi*Y(t-1) + e(t), e(t)~iidN(0,sig2)

clear;
clc;
addpath('C:\Users\이건희\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: DGP
T = 1000; 
c = 1;
phi = 0.8; 
sig2 = 1; 

uncondm = c/(1-phi);
uncondvar = sig2/(1-phi^2);

Ym = zeros(T,1); 
Ym(1) = uncondm + sqrt(uncondvar)*randn(1,1); 

for t = 2:T
    Ym(t) = c + phi*Ym(t-1) + sqrt(sig2)*randn(1,1);
end 

%% Step 2: Estimation 
Y = Ym(2:end);             % T by 1
YL = Ym(1:end-1);
T = rows(Y);
X = [ones(T,1) YL];      % T by 2
k = cols(X);

% Estimation
% parameter (c, phi, sig2)
theta0 = [0;0.5;0.1];
Data = [Y X];

% index
index = [1;2;3];
printi = 1;

% Optimization
[thetamx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);

% Compute t-value and p-value
diag_cov = diag(Vj);
stde = sqrt(diag_cov);
t_val = thetamx./stde; 
p_val = 2*(1 - cdf('t', abs(t_val), T-k)); 

%% Step 3: Results %%
disp('=========================================');
disp(['  Estimates ','  t value ', '  p value ']); 
disp('=========================================');
disp([thetamx t_val p_val]); 
disp('=========================================');

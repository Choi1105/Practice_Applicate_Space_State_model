%% Maximum likelihood estimation
% Non Linear Model
% y(t) = (1-rho)*ibar*x1(t) + (1-rho)*alpha*x2(t) + (1-rho)*beta*x3(t)...
%       + rho*x4(t) + e(t), e(t)~iidN(0,sig2)

clear;
clc;
addpath('C:\Users\Kim Young Min\Dropbox\Code\MATLAB\Lectures\M_library');

%% Step 1: DGP %%
T = 1000;

rho = 0.5; 
ibar = 4;
alpha = 2;
beta = -1; 
sig2 = 0.5; 

x1m = ones(T,1);
x2m = rand(T,1)*5;
x3m = rand(T,1)*5;
x4m = rand(T,1)*5;
em = sqrt(sig2)*randn(T,1); 

ym = (1-rho)*ibar*x1m + (1-rho)*alpha*x2m + (1-rho)*beta*x3m + rho*x4m + em;

%% Step 2: Estimation %%
% Data
Y = ym; 
X = [x1m x2m x3m x4m]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
theta0 = [0.1;0;0;0;0.1];

% index
index = [1;2;3;4;5];
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

          
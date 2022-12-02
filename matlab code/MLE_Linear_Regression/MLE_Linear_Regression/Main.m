%% Maximum likelihood estimation
% Linear Model
% Y = b1 + x2*b2 + e, Y=Xb+e

clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: DGP %%
T = 500;

b1 = 1; 
b2 = 2; 
b3 = 2.5;

sig2 = 2; 

x1m = ones(T,1);
x2m = rand(T,1)*5;
x3m = rand(T,1)*5;

em = sqrt(sig2)*randn(T,1); 

ym = x1m*b1 + x2m*b2 + x3m*b3 + em; 

%% Step 2: Estimation %%
% Data
Y = ym; 
X = [x1m x2m x3m]; 
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
psi0 = [0;0;0;1];

% index
index = [1;2;3;4];
printi = 1;

% Optimization
[psimx, fmax, Vj, Vinv] = SA_Newton(@lnlik,@paramconst,psi0,Data,printi,index);

% Compute t-value and p-value
diag_cov = diag(Vj);
stde = sqrt(diag_cov);
t_val = psimx./stde; 
p_val = 2*(1 - cdf('t', abs(t_val), T-k)); 

%% Step 3: Results %%
disp('=========================================');
disp(['  Estimates ','  t value ', '  p value ']); 
disp('=========================================');
disp([psimx t_val p_val]); 
disp('=========================================');

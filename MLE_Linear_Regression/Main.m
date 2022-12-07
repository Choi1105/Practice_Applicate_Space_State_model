%% Maximum likelihood estimation
% Linear Model
% Y = b1 + x2*b2 + e, Y=Xb+e

clear;
clc;
addpath('F:\Dropbox\Code\MATLAB\Lectures\M_library');

%% Step 1: DGP %%
T = 500;

b1 = 1; 
b2 = 2; 
sig2 = 3; 

x1m = ones(T,1);
x2m = rand(T,1)*5;
em = sqrt(sig2)*randn(T,1); 

ym = x1m*b1 + x2m*b2 + em; 

%% Step 2: Estimation %%
% Data
Y = ym; 
X = [x1m x2m];
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
psi0 = [0;0;1]; 

% index
index = [1;2;3];
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

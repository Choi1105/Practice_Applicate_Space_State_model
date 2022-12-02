%% Maximum likelihood estimation
% Autocorrelation
% Y(t) = b1 + b2*x2(t) + e(t)
% e(t) = rho*e(t-1) + v(t), v(t)


clear;
clc;
addpath('C:\Users\이건희\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: DGP %%
T = 500;
sig2 = 1;
rho = 0.5;
b1 = 0.3;
b2 = 0.2;

x1m = ones(T,1);
x2m = randn(T,1)*2;

em = zeros(T,1);
ym = zeros(T,1);
eL = sqrt(sig2/(1-rho))*randn(1,1);
for t = 1:T
   
    et = rho*eL + sqrt(sig2)*randn(1,1);     
    yt = b1*x1m(t) + b2*x2m(t) + et;
    
    em(t) = et;
    ym(t) = yt;
    
    eL = et;
end


%% Step 2: Estimation %%
% Data
Y = ym; 
X = [x1m x2m]; 

Data = [Y X];

T = rows(X);
k = cols(X);

% initial value 
theta0 = [0;0;0.1;0.1];

% index
index = [1;2;3;4];
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

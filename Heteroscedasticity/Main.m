%% Maximum likelihood estimation
% Heteroscedasticity
% Y(t) = b1 + b2*x2(t) + e(t), e(t)~iidN(0,sig2t)
% sig2t = a0 + a1*z(t)

clear;
clc;
addpath('C:\Users\Kim Young Min\Dropbox\Code\MATLAB\Lectures\M_library');

%% Step 1: DGP %%
T = 1000;

b1 = 1; 
b2 = 2;

a0 = 0.5;
a1 = 5;

x1m = ones(T,1);
x2m = randn(T,1)*2;
zm = rand(T,1)*5;

em = zeros(T,1);
ym = zeros(T,1);

for t = 1:T
    sig2t = a0 + a1*zm(t);
    et = sqrt(sig2t)*randn(1,1);     
    yt = b1*x1m(t) + b2*x2m(t) + et;
    
    em(t) = et;
    ym(t) = yt;
end


%% Step 2: Estimation %%
% Data
Y = ym; 
X = [x1m x2m]; 
Z = zm;
Data = [Y X Z];

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

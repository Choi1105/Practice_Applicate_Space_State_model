%% Maximum likelihood estimation
% Non Linear Model
% y(t) = (1-rho)*ibar*x1(t) + (1-rho)*alpha*x2(t) + (1-rho)*beta*x3(t)...
%       + rho*x4(t) + e(t), e(t)~iidN(0,sig2)

clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% Step 1: DGP %%
ym = readmatrix('KOSPI','sheet','데이터','range','B2:B110');              % KOSPI 지수 변화율
x1m = ones(109,1);
x2m = readmatrix('KOSPI','sheet','데이터','range','C2:C110');             % PMI 변화율
x3m = readmatrix('KOSPI','sheet','데이터','range','D2:D110');             % 환율 변화율
x4m = readmatrix('KOSPI','sheet','데이터','range','F2:F110');             % 한국 10년 국채 금리

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

          
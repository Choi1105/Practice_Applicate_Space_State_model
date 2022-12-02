%% Maximum likelihood estimation
% Autocorrelation
% Y(t) = b1 + b2*x2(t) + e(t)
% e(t) = rho*e(t-1) + v(t), v(t)

clear;
clc;

%% Step 1: DGP %%
ym = readmatrix('한전데이터','sheet','Sheet1','range','I4:I67');   % 한전 영업이익 차분

x1m = ones(64,1);
x2m = readmatrix('한전데이터','sheet','Sheet1','range','M4:M67');  % LNG 정산단가

rho = readmatrix('한전데이터','sheet','Sheet1','range','P4:P67');  % 신재생 에너지 발전량


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

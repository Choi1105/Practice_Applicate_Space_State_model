%% Maximum likelihood estimation 
% Hypothesis, Likelihood Ratio Test
% Unrestricted Model
% Y(t) = b1 + b2*x2(t) + e(t), e(t)~iidN(0,sig2)

% Restricted Model
% Y(t) = b1 + e(t), e(t)~iidN(0,sig2)


clear;
clc;
%% Step 1: DGP %%
T = 1000;

b1 = 1; 
b2 = 0; 
sig2 = 2; 

x1m = ones(T,1);
x2m = rand(T,1)*5;
em = sqrt(sig2)*randn(T,1); 

ym = x1m*b1 + x2m*b2 + em; 

%% Step 2: Unrestricted Regression
% Data
Y = ym; 
X = [x1m x2m]; 
Data = [Y X];

% initial value 
theta0 = [0;0;1];

% index
index = [1;2;3];
printi = 1;

% Optimization
[~, lnL_U, ~, ~] = SA_Newton(@lnlik,@paramconst,theta0,Data,printi,index);

%% Step 3: Restricted Regression
X = x1m; 
Data = [Y X];

% initial value 
theta0 = [0;1];

% index
index = [1;2];
printi = 1;

% Optimization
[~, lnL_R, ~, ~] = SA_Newton(@lnlik_R,@paramconst_R,theta0,Data,printi,index);

%% Step 4: LR test
nR = 1;
LR = -2*(lnL_R - lnL_U);
p_val = 1 - cdf('Chi',LR,nR);

disp('================================================='); 
disp(' LR test for X2'); 
disp('-------------------------------------------------'); 
disp([' H0      ', 'Test stat.', '   P_Value']);
disp('-------------------------------------------------'); 
disp([' X2      ', num2str(LR), '      ', num2str(p_val)]); 
disp('-------------------------------------------------'); 

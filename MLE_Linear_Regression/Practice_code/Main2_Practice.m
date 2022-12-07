%% Maximum likelihood estimation
% Linear Model
% Y = b1 + x2*b2 + e, Y=Xb+e
% Practice Question : Add the variable "x3"
% Y = b1 + x2*b2 + x3*b3 + e, Y=Xb+e

clear;
clc;
addpath('F:\Dropbox\Code\MATLAB\Lectures\M_library');

%% Step 1: DGP %%
T = 500;

b1 = 1; 
b2 = 2; 
b3 = 3;
sig2 = 2; 

x1m = ones(T,1);
x2m = rand(T,1)*5;
x3m = rand(T,1)*5; % Dimension (T,1) (To add variable X3
em = sqrt(sig2)*randn(T,1); 

ym = x1m*b1 + x2m*b2 + x3m*b3 + em; 

%% Step 2: Estimation %%
% Data
Y = ym; 
X = [x1m x2m x3m]; % Variable add in this line
Data = [Y X];
T = rows(X);
k = cols(X);

% initial value 
psi0 = [0;0;0;1]; 
% psi0's element means that [Mu(beta1) , Beta2 , Beta3(add variable) , sig2]
% in this line you must revise the paramconst.m file (line num 6) theta("n")
% "n" means that sig2 that "n"th element psi0 vector. 
% so, In this case you revise the "n" is 4
% And, the same reason, you revise the lnlik.m file's line(11~13) too 

% index
index = [1;2;3;4]; 
% the same reason, you revise this line too (Generate the New index)

printi = 1; 
% print option 0, 1, 2 -> Try it (maybe the option exist more than 3?..Well..I don't know)

% Optimization
[psimx, fmax, Vj, Vinv] = SA_Newton(@lnlik2_Practice,@paramconst2_Practice,psi0,Data,printi,index);
% At this code you set the wrong number of "n", index or something that I
% refer in before Error will reported So, If the Error reported that after
% run the code you have to revise something codes.


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

% How can I read the result?

% Estimates is M, beta2, beta3, sig2 
% If you revise the code in right way 
% you can get the similar result that you set at line 14~17  

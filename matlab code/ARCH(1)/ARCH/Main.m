%% Heteroskedasticity Model
% ARCH 

% Model Specification
% y(t) = mu + e(t), e(t) = sig(t)*v(t), v(t)~iidN(0,1)
% sig(t)^2 = alpha0 + alpha1*e(t-1)^2

clc; 
clear;

%% Step1 :Data Generating Process
% sample size
T = 1000; 

% True parameter
mu = 0.5;
alpha0 = 0.3;
alpha1 = 0.8;
tru_para = [mu;alpha0;alpha1];

% Pre-allocation
ym = zeros(T,1);
em = zeros(T,1);
tru_var = zeros(T,1);

% Unconditional Initial Volatility
tru_var(1) = alpha0/(1-alpha1); 

% Initial error/dependent variable
em(1) = sqrt(tru_var(1))*randn(1,1);
ym(1) = mu + em(1);

% Generate order: Vol(t),  e(t), y(t)
for t = 2:T

    tru_var(t) = alpha0 + alpha1* (em(t-1)^2); 
    em(t) = sqrt(tru_var(t))*randn(1,1); 
    ym(t) = mu + em(t);

end

%% Step 2: Maxmimum Likelihood Estimation
% Data
Data = ym;

% Initial values
psi0 = [0;0.1;0.2];

% Index
indbj = [1;2;3];

% printi = 1 => See the opimization produdure
% printi = 0 => NOT see the opimization produdure
printi = 1; 

% Optimization
[thetamx, fmax, V, Vinv] = SA_Newton(@lnlik,@paramconst,psi0,Data,printi,indbj);

diag_V = diag(V);
stde = sqrt(diag_V);
t_val = thetamx./stde;
p_val = 2*(1 - cdf('t',abs(t_val),T-1));

%% Step 3: Table / Figure Results 
% Table
disp('===========================================================');
disp(['    Index ','  True Para ', ' Estimates ', ' t value ',  ' p value']);
disp('===========================================================');
disp([indbj tru_para thetamx t_val p_val]);
disp('===========================================================');

% Figure
% Compute esitamted variance and residual
[~, Varm , Residm, Std_Residm] = lnlik(thetamx,Data);

% Index 
i = 1:rows(Varm); 
i = i';

figure
plot(i, tru_var ,'k', i, Varm, 'b:', 'LineWidth',1.5);
legend('True','Filtered'); 
title('True and Estimated');

figure
plot(i, Residm ,'r', i, Std_Residm, 'b:', 'LineWidth',1.5);
legend('Resid','Std Resid'); 
title('Residual and Standard Residual');

% Test Engle's ARCH LM test
[phi_hat, sig2_hat, R2] = OLS_ARp(Residm.^2, 1);
DoF = 1;
LM = (T-1)*R2;
p_val = 1 - cdf('Chi',LM,DoF);
disp('================================================='); 
disp(' ARCH LM test'); 
disp('-------------------------------------------------'); 
disp([' H0         ', '   Test stat.', '   P_Value']);
disp('-------------------------------------------------'); 
disp([' No Hetero      ', num2str(LM), '      ', num2str(p_val)]); 
disp('-------------------------------------------------'); 


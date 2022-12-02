%% Heteroskedasticity Model
% GARCH in Mean (GARCHM)

% Model Specification
% y(t) = mu + beta*sig(t) + e(t), e(t)|I(t-1) = sig(t)*v(t), v(t)~iid N (0,1)
% sig(t)^2 = alpha0 + alpha1*e(t-1)^2 + gam1*sig(t-1)^2
clear; 
clc;   

%% Step1 :Data Generating Process
% sample size
T = 1000;

% true value
b0 = 2;
b1 = 1;
a0 = 1;
a1 = 0.5;
gam1 = 0.3;
tru_para = [b0;b1;a0;a1;gam1];

% Pre-allocation
ym = zeros(T,1); 
em = zeros(T,1); 
tru_var = zeros(T,1); 

% Unconditional Initial Volatility
tru_var(1) = a0/(1-a1-gam1); 

% Initial error/dependent variable
em(1) = randn(1,1)*sqrt(tru_var(1)); 
ym(1) = b0 + b1*sqrt(tru_var(1)) + em(1); 

% Generate order: Vol(t),  e(t), y(t)
for t = 2:T

    tru_var(t) = a0 + a1*(em(t-1)^2) + gam1*tru_var(t-1); 
    em(t) = sqrt(tru_var(t))*randn(1,1); 
    ym(t) = b0 + b1*sqrt(tru_var(t)) + em(t);

end


%% Step 2: Maxmimum Likelihood Estimation
% Data
Data = ym;

% Initial values
psi0 = [b0;b1;a0;a1;gam1];

% Index
indbj = [1;2;3;4;5];

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
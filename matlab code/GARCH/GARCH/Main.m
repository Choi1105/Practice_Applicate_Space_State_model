%% Heteroskedasticity Model
% GARCH
% Estimated by MLE

% Model Specification
% y(t) = mu + e(t), e(t)|I(t-1) = sig(t)*v(t), v(t)~iid N (0,1)
% sig(t)^2 = alpha0 + alpha1*e(t-1)^2 + gam1*sig(t-1)^2

clear; 
clc;   

%% Step1 :Data Generating Process
% sample size
T = 1000; 

% True parameter
b1 = 0.1;
b2 = 0.1;
b3 = -0.2;
b4 = 0.2;
a0 = 0.2;
a1 = 0.5;
gam1 = 0.3;
tru_para = [b1;b2;b3;b4;a0;a1;gam1];

x1m = ones(T,1);
x2m = 5*rand(T,1);
x3m = 5*rand(T,1);
x4m = 5*rand(T,1);

% Pre-allocation
ym = zeros(T,1); 
em = zeros(T,1); 
tru_var = zeros(T,1); 

% Unconditional Initial Volatility
tru_var(1) = a0/(1-a1-gam1); 

% Initial error/dependent variable
em(1) = randn(1,1)*sqrt(tru_var(1)); 
ym(1) = x1m(1)*b1 + x2m(1)*b2 + x3m(1)*b3 + x4m(1)*b4 + em(1); 

% Generate order: Vol(t),  e(t), y(t)
for t = 2:T 

    tru_var(t) = a0 + a1* (em(t-1)^2) + gam1*tru_var(t-1); 
    em(t) = sqrt(tru_var(t))*randn(1,1); 
    ym(t) = x1m(t)*b1 + x2m(t)*b2 + x3m(t)*b3 + x4m(t)*b4 + em(t) ;

end

%% Step 2: Maxmimum Likelihood Estimation
% Data
Data = [ym, x1m, x2m, x3m, x4m];

% Initial values
psi0 = [0;0;0;0;0.1;0.1;0.1];

% Index
indbj = [1;2;3;4;5;6;7];

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

sig2_Pred = thetamx(5) + thetamx(6)*Residm(T) + thetamx(7)*Varm(T)
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
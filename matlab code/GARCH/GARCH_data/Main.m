%% Heteroskedasticity Model
% GARCH
% Estimated by MLE

% Model Specification
% y(t) = X(t)'*beta + e(t), e(t)|I(t-1) = sig(t)*v(t), v(t)~iid N (0,1)
% sig(t)^2 = alpha0 + alpha1*e(t-1)^2 + gam1*sig(t-1)^2

clear; 
clc;   

%% Step1 :Loading Data
% Period : 2003년 1월 ~ 2022년 5월
ym = xlsread('Kospi.csv','Kospi','C4:C5560');

T = rows(ym);


x1m = ones(T,1);
x2m = xlsread('Kospi.csv','Kospi','C3:C5559');   

%% Step 2: Maxmimum Likelihood Estimation
% Data
Data = [ym, x1m, x2m];

% Initial values
psi0 = [0;0;0.1;0.1;0.1];

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
disp(['    Index ', ' Estimates ', ' t value ',  ' p value']);
disp('===========================================================');
disp([indbj thetamx t_val p_val]);
disp('===========================================================');

% Figure
% Compute esitamted variance and residual
[~, Varm , Residm, Std_Residm] = lnlik(thetamx,Data);

sig2_Pred = thetamx(3) + thetamx(4)*Residm(T) + thetamx(5)*Varm(T);

% Index 
i = 1:rows(Varm); 
i = i';

figure
plot(i, ym, 'k', i, Varm, 'b:', 'LineWidth',1.5);
legend('Data', 'Filtered Variance'); 
title('Data and Estimated Variance');


figure
plot(i, Residm ,'r', i, Std_Residm, 'b:', 'LineWidth',1.5);
legend('Resid','Std Resid'); 
title('Residual and Standard Residual');
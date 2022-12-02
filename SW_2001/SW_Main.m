%% Recursive VAR model 
% Using Cholesky Decomposition
% Replicate Stock & Watson(2001) VAR(4)
% Shortrun restriction

clear;
clc;

p = 4;  % lag
J = 24; % impulse response 

%% Data Loading
% Inf, Un, FFR
% Sample period: 1960Q1 ~ 2001Q3

[Data, ~, ~] = xlsread('watson.xlsx', 'var', 'A1:C167');

T = rows(Data); % Sample size

inf1 = Data(:,1); % Inflation Rate
un1 = Data(:,2);  % Unemployement
ffr1 = Data(:,3); % Federal Funds Rate

% Deviation from mean
inf = inf1 - mean(inf1);
un = un1 -  mean(un1);
ffr = ffr1 -  mean(ffr1);

% Data
y = [inf un ffr]; % T by k
k = cols(y);

% Y, (T-p) by k
y0 = y(p+1:T,:); 

% Lag
y1 = y(p:T-1,:);  
y2 = y(p-1:T-2,:);
y3 = y(p-2:T-3,:);
y4 = y(p-3:T-4,:);

% X
y_L = [y1 y2 y3 y4]; 

% 
T = rows(y0);

% inv(X'X)X'Y(OLS), p*k by k
phihat = (y_L'*y_L)\(y_L'*y0); 

% Making F matrix,  p*k by p*k
F = [phihat';eye((p-1)*k), zeros((p-1)*k,k)];

% Making Omega matrix
uhat = y0 - y_L*phihat; % T-p by k
Omegahat = uhat'*uhat/(T-p*k); % k by k

%% Estimation of inv(B) Using Cholesky Decomposition
invB = chol(Omegahat)';  % Lower triangular matrix

%% Impulse Response 
FF = eye(p*k);

% Impulse response 
theta_inf = zeros(J+1,k); % Impulse response of INF to each shock 
theta_un = zeros(J+1,k);  % Impulse response of UN to each shock 
theta_ffr = zeros(J+1,k); % Impulse response of FFR to each shock 

for j = 1:J+1 % j=1 은 사실 j=0일때 이다.
    
    psi_j = FF(1:k, 1:k);       % First black of F^j Matrix, k by k
    theta  = psi_j*invB;        % Impulse Response Coef, k by k 
    
    theta_inf(j,:) = theta(1,:); % 1 by k 각 Shock에 대한 Inf 반응
    theta_un(j,:) = theta(2,:);  % 1 by k 각 Shock에 대한 Un 반응 
    theta_ffr(j,:) = theta(3,:); % 1 by k 각 Shock에 대한 FFR 반응 
    
    FF = FF*F;                   % F->F^2->F^3...loop
end
   
%% Impulse Response Plot
zeroline = zeros(J+1,1);
figure
subplot(3,3,1);
plot([theta_inf(:,1) zeroline]);
title('Inflation shock to Inflation');
subplot(3,3,2);
plot([theta_un(:,1) zeroline]);
title('Inflation shock to Un');
subplot(3,3,3);
plot([theta_ffr(:,1) zeroline]);
title('Inflation shock to FFR');

subplot(3,3,4);
plot([theta_inf(:,2) zeroline]);
title('Un shock to Inflation');
subplot(3,3,5);
plot([theta_un(:,2) zeroline]);
title('Un shock to Un');
subplot(3,3,6);
plot([theta_ffr(:,2) zeroline]);
title('Un shock to FFR');

subplot(3,3,7);
plot([theta_inf(:,3) zeroline]);
title('FFR shock to Inflation');
subplot(3,3,8);
plot([theta_un(:,3) zeroline]);
title('FFR shock to Un');
subplot(3,3,9);
plot([theta_ffr(:,3) zeroline]);
title('FFR shock to FFR');

%% Variance decomposition
vd_inf = zeros(J+1,k); 
vd_un = zeros(J+1,k);
vd_ffr = zeros(J+1,k);

theta_inf2 = theta_inf.*theta_inf;
theta_un2 = theta_un.*theta_un;
theta_ffr2 = theta_ffr.*theta_ffr;

for i = 1:k

    var_inf = 0;
    var_un = 0;
    var_ffr = 0;

    for j = 1:J+1
    
        % V.D. of inflation rate
        var_inf = var_inf + sum(theta_inf2(j,:));
        vd_inf(j,i) = sum(theta_inf2(1:j,i))/var_inf;

        % V.D. of Un
        var_un = var_un + sum(theta_un2(j,:));
        vd_un(j,i) = sum(theta_un2(1:j,i))/var_un;

        % V.D. of FFR   
        var_ffr = var_ffr + sum(theta_ffr2(j,:));
        vd_ffr(j,i) = sum(theta_ffr2(1:j,i))/var_ffr;

    end
      
end

% Variance Decomposition Plot
figure
zeroline = zeros(J+1,1);
subplot(3,3,1);
plot([vd_inf(:,1) zeroline ones(J+1,1)])
title('Inflation shock to Inflation')
subplot(3,3,2);
plot([vd_inf(:,2) zeroline ones(J+1,1)])
title('Un shock to inflation')
subplot(3,3,3);
plot([vd_inf(:,3) zeroline ones(J+1,1)])
title('FFR shock to inflation')

subplot(3,3,4);
plot([vd_un(:,1) zeroline ones(J+1,1)])
title('Inflation shock to Un')
subplot(3,3,5);
plot([vd_un(:,2) zeroline ones(J+1,1)])
title('Un shock to Un')
subplot(3,3,6);
plot([vd_un(:,3) zeroline ones(J+1,1)])
title('FFR shock to Un')

subplot(3,3,7);
plot([vd_ffr(:,1) zeroline ones(J+1,1)])
title('Inflation shock to FFR')
subplot(3,3,8);
plot([vd_ffr(:,2) zeroline ones(J+1,1)])
title('Un shock to FFR')
subplot(3,3,9);
plot([vd_ffr(:,3) zeroline ones(J+1,1)])
title('FFR shock to FFR')

% Model Selection
AIC = log(det(Omegahat)) + 2*k^2*p/T;
BIC = log(det(Omegahat)) + k^2*p*log(T)/T;

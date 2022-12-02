%% Recursive VAR model 
% Using Cholesky Decomposition
% Replicate Stock & Watson(2001) VAR(4)
% Shortrun restriction

clear;
clc;

p = 2;  % lag
J = 20; % impulse response 

%% Data Loading
% Inf, Un, FFR
% Sample period: 1960Q1 ~ 2001Q3

[Data, ~, ~] = xlsread('VAR_data.xlsx', 'Sheet3', 'A1:D1238');

T = rows(Data); % Sample size

KOS1 = Data(:,1);  % KOSPI 종합
Sm1 = Data(:,2);   % KOSPI 소형주
Mm1 = Data(:,3);   % KOSPI 중형중
Bm1 = Data(:,4);   % KOSPI 대형주



% Deviation from mean
KOS = KOS1 - mean(KOS1);
Sm = Sm1 -  mean(Sm1);
Mm = Mm1 -  mean(Mm1);
Bm = Bm1 - mean(Bm1);



% Data
y = [KOS Sm Mm Bm]; % T by k
k = cols(y);

% Y, (T-p) by k
y0 = y(p+1:T,:); 

% Lag
y1 = y(p:T-1,:);  
y2 = y(p-1:T-2,:);

% X
y_L = [y1 y2]; 

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
theta_KOS = zeros(J+1,k); % Impulse response of KOSPI to each shock 
theta_Sm = zeros(J+1,k);  % Impulse response of Small Stock to each shock 
theta_Mm = zeros(J+1,k);  % Impulse response of Medium Stock to each shock 
theta_Bm = zeros(J+1,k);


for j = 1:J+1 % j=1 은 사실 j=0일때 이다.
    
    psi_j = FF(1:k, 1:k);       % First black of F^j Matrix, k by k
    theta  = psi_j*invB;        % Impulse Response Coef, k by k 
    
    theta_KOS(j,:) = theta(1,:); % 1 by k 각 Shock에 대한 Inf 반응
    theta_Sm(j,:) = theta(2,:);  % 1 by k 각 Shock에 대한 Un 반응 
    theta_Mm(j,:) = theta(3,:); % 1 by k 각 Shock에 대한 FFR 반응 
    theta_Bm(j,:) = theta(4,:);

    FF = FF*F;                   % F->F^2->F^3...loop
end
   
%% Impulse Response Plot
zeroline = zeros(J+1,1);
figure
subplot(4,4,1);
plot([theta_KOS(:,1) zeroline]);
title('Kos shock to Kos');
subplot(4,4,2);
plot([theta_Sm(:,1) zeroline]);
title('Kos shock to Sm');
subplot(4,4,3);
plot([theta_Mm(:,1) zeroline]);
title('Kos shock to Me');
subplot(4,4,4);
plot([theta_Bm(:,1) zeroline]);
title('Kos shock to Bm');


subplot(4,4,5);
plot([theta_KOS(:,2) zeroline]);
title('Sm shock to Kos');
subplot(4,4,6);
plot([theta_Sm(:,2) zeroline]);
title('Sm shock to Sm');
subplot(4,4,7);
plot([theta_Mm(:,2) zeroline]);
title('Sm shock to Mm');
subplot(4,4,8);



subplot(4,4,9);
plot([theta_KOS(:,3) zeroline]);
title('Mm shock to Kos');
subplot(4,4,10);
plot([theta_Sm(:,3) zeroline]);
title('Mm shock to Sm');
subplot(4,4,11);
plot([theta_Mm(:,3) zeroline]);
title('Mm shock to Mm');
subplot(4,4,12);
plot([theta_Bm(:,3) zeroline]);
title('Mm shock to Bm');

subplot(4,4,13);
plot([theta_KOS(:,4) zeroline]);
title('Bm shock to Kos');
subplot(4,4,14);
plot([theta_Sm(:,4) zeroline]);
title('Bm shock to Sm');
subplot(4,4,15);
plot([theta_Mm(:,4) zeroline]);
title('Bm shock to Mm');
subplot(4,4,16);
plot([theta_Bm(:,4) zeroline]);
title('Bm shock to Bm');

%% Variance decomposition
vd_KOS = zeros(J+1,k); 
vd_Sm = zeros(J+1,k);
vd_Mm = zeros(J+1,k);
vd_Bm = zeros(J+1,k);

theta_KOS2 = theta_KOS.*theta_KOS;
theta_Sm2 = theta_Sm.*theta_Sm;
theta_Mm2 = theta_Mm.*theta_Mm;
theta_Bm2 = theta_Bm.*theta_Bm;


for i = 1:k

    var_KOS = 0;
    var_Sm = 0;
    var_Mm = 0;
    var_Bm = 0;

    for j = 1:J+1
    
        % V.D. of KOSPI
        var_KOS = var_KOS + sum(theta_KOS2(j,:));
        vd_KOS(j,i) = sum(theta_KOS2(1:j,i))/var_KOS;

        % V.D. of Sm
        var_Sm = var_Sm + sum(theta_Sm2(j,:));
        vd_Sm(j,i) = sum(theta_Sm2(1:j,i))/var_Sm;

        % V.D. of Mm   
        var_Mm = var_Mm + sum(theta_Mm2(j,:));
        vd_Mm(j,i) = sum(theta_Mm2(1:j,i))/var_Mm;

        % V.D. of Bm
        var_Bm = var_Bm + sum(theta_Bm2(j,:));
        vd_Bm(j,i) = sum(theta_Bm2(1:j,i))/var_Bm;


    end
      
end

% Variance Decomposition Plot
figure
zeroline = zeros(J+1,1);
subplot(4,4,1);
plot([vd_KOS(:,1) zeroline ones(J+1,1)])
title('Kos shock to Kos');
subplot(4,4,2);
plot([vd_KOS(:,2) zeroline ones(J+1,1)])
title('Kos shock to Sm');
subplot(4,4,3);
plot([vd_KOS(:,3) zeroline ones(J+1,1)])
title('Kos shock to Mm');
subplot(4,4,4);
plot([vd_KOS(:,4) zeroline ones(J+1,1)])
title('Kos shock to Bm');


subplot(4,4,5);
plot([vd_Sm(:,1) zeroline ones(J+1,1)])
title('Sm shock to KOS')
subplot(4,4,6);
plot([vd_Sm(:,2) zeroline ones(J+1,1)])
title('Sm shock to Sm')
subplot(4,4,7);
plot([vd_Sm(:,3) zeroline ones(J+1,1)])
title('Sm shock to Mm')
subplot(4,4,8);
plot([vd_Sm(:,4) zeroline ones(J+1,1)])
title('Sm shock to Bm')

subplot(4,4,9);
plot([vd_Mm(:,1) zeroline ones(J+1,1)])
title('Mm shock to KOS')
subplot(4,4,10);
plot([vd_Mm(:,2) zeroline ones(J+1,1)])
title('Mm shock to Sm')
subplot(4,4,11);
plot([vd_Mm(:,3) zeroline ones(J+1,1)])
title('Mm shock to Mm')
subplot(4,4,12);
plot([vd_Mm(:,4) zeroline ones(J+1,1)])
title('Mm shock to Bm')

subplot(4,4,13);
plot([vd_Bm(:,1) zeroline ones(J+1,1)])
title('Bm shock to KOS')
subplot(4,4,14);
plot([vd_Bm(:,2) zeroline ones(J+1,1)])
title('Bm shock to Sm')
subplot(4,4,15);
plot([vd_Bm(:,3) zeroline ones(J+1,1)])
title('Bm shock to Mm')
subplot(4,4,16);
plot([vd_Bm(:,4) zeroline ones(J+1,1)])
title('Bm shock to Bm')

% Model Selection
AIC = log(det(Omegahat)) + 2*k^2*p/T;
BIC = log(det(Omegahat)) + k^2*p*log(T)/T;

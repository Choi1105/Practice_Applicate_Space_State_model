%% Heteroskedasticity Model
% Stochastic Volatility
% Gibbs Sampler Based on Kim, Shephard, and Chib (1998, ReStud)
% Application: Oil price

% Model Specification
% y(t) = exp(h(t)/2)e(t), e(t)~iidN(0,1)
% h(t) = mu + phi*h(t-1) + u(t),u(t)~iidN(0,Sig2)

% Model Transformation State-space model with Normal Approximation
% Measurement Eq.
% y(t)_Star = h(t) + e(t)_Star, e(t)_Star|s(t) ~ N(m_s(t), v_s(t))
% where y(t)_Star = ln(y(t)^2), e(t)_Star= ln(e(t)^2)

% Transition Eq.
% h(t) = mu + phi*h(t-1) + u(t),u(t)~iidN(0,Sig2)

clear
clc

%% Step 1: Date Loading
% Data Loading
Y = xlsread('Oil price data', 'Oil price', 'K230:K408');    % Oil Price (t)

% Sample size
T = rows(Y);

x1m = ones(T,1);
x2m = xlsread('원유생산및소비', 'sheet1', 'C2:C180');        % Oil Production(t-1)
x3m = xlsread('Oil price data', 'Oil price', 'K229:K407');  % Oil Price (t-1)
x4m = xlsread('Oil price data', 'Oil price', 'M229:M407');  % Oil Consumption (World Industrial Price Index) (t-1)

X = [x1m x2m x3m x4m];

% Demean
ym = demeanc(Y);

%% Step 3: Estimation by MCMC Method: Gibbs Sampler
% Data transformation
ysm = log(ym.*ym);

% MCMC size 
n0 = 1000;   % burn-in period
n1 = 2000;   % MCMC size after burn-in period
freq = 100;  % frequeccy of reporting intermediate
n = n0 + n1; % Sampling in total

% Prior distribution
% Normal for mu and phi
mu_phi_M = 0.1*ones(2,1);
mu_phi_V = 1*eye(2);

% precision matrix
precb_ = inv(mu_phi_V);

% Inverse Gamma for sig2
s2m_ = 1;  
s2se_ = 2;  
[v0,d0] = paramig(s2m_,s2se_);

% Selection of the mixing distribution to be log(Chi-square(1))
% Refer to Kim, Shephard, and Chib (1998, ReStud), Table 4
pm = zeros(7,1);
pm(1) = 0.00730;
pm(2) = 0.10556;
pm(3) = 0.00002;
pm(4) = 0.04395;
pm(5) = 0.34001;
pm(6) = 0.24566;
pm(7) = 0.25750;

msm = zeros(7,1);
msm(1) = -11.40039;
msm(2) = -5.24321;
msm(3) = -9.83726;
msm(4) = 1.50746;
msm(5) = -0.65098;
msm(6) = 0.52478;
msm(7) = -2.35859;

vsm = zeros(7,1);
vsm(1) = 5.79596;
vsm(2) = 2.61369;
vsm(3) = 5.17950;
vsm(4) = 0.16735;
vsm(5) = 0.64009;
vsm(6) = 0.34023;
vsm(7) = 1.26261;

% Initial values
mu_phi = mu_phi_M;
sig2_inv = 1/s2m_;
hm = ysm;

% Pre-allocation
Hm = zeros(n,T);
Volm = zeros(n,T);
Sig2m = zeros(n,1);
mu_phim = zeros(n,2);

for iter = 1:n
% Step 0: Data
    Y = hm(2:end); % h(t)
    X = [ones(rows(hm)-1,1) hm(1:end-1)]; % [1, h(t-1)]
    
% Step 1: Mu and Phi sammpling
    [mu_phi] = Gen_mu_phi(Y,X,mu_phi_M,precb_,sig2_inv,mu_phi);
    mu_phim(iter,:) = mu_phi'; 

% Step 2: Sig2 sampling
    [sig2_inv,sig2] = Gen_Sigma(Y,X,v0,d0,mu_phi);
    Sig2m(iter,:) = sig2;

% Step 3: Sm Sampling 
    sm = Gen_Sm(ysm,hm,pm,msm,vsm);

% Step 4: Hm Sampling
    [hm,~] = Gen_Fm(ysm,mu_phi,sig2,sm,msm,vsm);     
    Hm(iter,:) = hm';  
    
% Compute volatility
    vol = exp(hm/2); 
    Volm(iter,:) = vol';
    
% Report intermediate result     
    if isequal(floor(iter/freq),iter/freq) == 1 
        prt2(mu_phi,sig2,iter);      
    end    
       
end
    

%% Step 4: Results
% Burn-in
MHm = [mu_phim Sig2m];
MHm = MHm(n0+1:n,:);
Volm = Volm(n0+1:n,:);
Hm = Hm(n0+1:n,:);

% Summary of posterior distribution
alpha = 0.025;    % credible interval
maxac = 200;      % inefficiency factor lag
is_postplot = 1;  % Plot posterior or NOT

postmom = MHout(MHm,alpha,maxac,is_postplot); 

% Table 
disp('MCMC Result');
disp('-----------------------------------------------------------');
disp(['     Mean ', '     2.5% ', '     97.5% ',  '   ineff.   ',' p-value']);
disp('-----------------------------------------------------------');
disp([postmom(:,2) postmom(:,4) postmom(:,6) postmom(:,7) postmom(:,8)]);
disp('===========================================================');

% Figure
% Posterior mean of Volatility
volm = meanc(Volm);

% Index 
i = 1:rows(volm); 
i = i';

% 2008/07/07 ~ 2014/07/03
figure
plot(i, ym ,'b:', i, volm, 'k', 'LineWidth',1.5);
legend('Oil price','Posterior mean'); 
title('Oil price and Estimated Volatility');
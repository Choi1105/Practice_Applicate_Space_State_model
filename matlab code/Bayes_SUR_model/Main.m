%% Seemingly Unrelated Regression Model
% Model 
% y(1t) = x(1t)'*beta1,St + e(1t)
% y(2t) = x(2t)'*beta2,St + e(2t)
% y(3t) = x(3t)'*beta3,St + e(3t)
% where x(1t) = x(2t) = x(3t) = ki by 1
%       (e(1t) e(2t) e(3t))' = e(t), where e(t) ~ N(0, Omega), p by 1
clear
clc

%% STEP 1: Data Generating Process
% Data information
% sample size(T), number of regression(p), 
% number of coefficient in each regression i(ki), total coefficient(k)
T = 1000; 

p = 3;
ki = 2;
k = p*ki;

% True beta
b1 = [1;2]; 
b2 = [3;4]; 
b3 = [5;6]; 
beta = [b1;b2;b3]; 

% True Omega = V*Gam*V, (p by p)
% where Gam = correlation Matrix / V = Volatility matrix (diagonal)
Gam = [1, 0.8, 0.8;0.8, 1, 0.8;0.8, 0.8, 1]; 
V = 0.5*eye(p);
Omega = V*Gam*V;  

% Exogenous variables in each regression
X1 = [ones(T,1) rand(T,1)*10];
X2 = [ones(T,1) rand(T,1)*10];
X3 = [ones(T,1) rand(T,1)*10];

% Pre-allocation for ym
ym = zeros(T,p);

for t = 1:T
    
    % at each time t, exogenous variable matrix, p by k
    xt = zeros(p,k);  
    xt(1,1:ki) = X1(t,:);
    xt(2,(ki+1):2*ki) = X2(t,:);
    xt(3,(2*ki+1):3*ki) = X3(t,:);
  
    % at each time t, dependent variables, p by 1 
    ym(t,:) = xt*beta + chol(Omega)'*randn(p,1);  
    
end

%% STEP 2: Set-up
% Data
Y = ym;
X = [X1 X2 X3];

% Burn in(n0) / Actual sampling size(n1)
n0 = 100; 
n1 = 1000; 

% Conjugate Prior Distributions
% Normal prior for beta(=vec(Phi))
b0 = zeros(k,1); 
var0 = 1*eye(k);

% Inverse Wishart prior for inv(Omega)
% Wishart dist, X ~ W(p)(R,v), E(X) = R*v
% Inverse Wishart dist, inv(X) ~ IW(p)(inv(R),v), E(X) = inv(R)/(v-k-1)
nu0 = 10;
R0 = invpd(eye(p))/(nu0-p-1);  

% Structure Variables
Spec.Y = Y;
Spec.X = X;

Spec.b0 = b0;
Spec.var0 = var0;
Spec.nu0 = nu0;
Spec.R0 = R0;

%% STEP 3: Gibbs Sampling
MHm = Gibbs_SUR(n0,n1,Spec);

%% STEP 4: Report Posterior Distributions
% criteria for output
alpha = 0.025;      % credible percent/2
maxac = 200;        % maximum number of autocorrelation
is_postplot = 1;    % = 1, display posterior distribution, = 0, NOT

% Posterior Moments / Sampling efficiency
postmom = MHout(MHm,alpha,maxac,is_postplot);



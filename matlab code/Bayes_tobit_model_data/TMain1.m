clc
clear
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library');

%% DGP
Data = xlsread('0519.xlsx','2019','G2:AW2289');

income = Data(:,1);
assets = Data(:,2);
wage = Data(:,4);
age = Data(:,7);
d_digital = Data(:,15);
benke = Data(:,11);
%% yuce
Y = benke;
T = rows(Y);
X = [ones(T,1) d_digital wage income assets age];
k = cols(X);


beta_0 = zeros(k,1);
b0 = eye(k);
b0inv = inv(b0);

v_0 = 10;
d_0 = 10;

n1 = 1000;
n0 = 1000;
n = n0 + n1;

Z = mean(Y); %Z 初始值
betam = zeros(n,k);
sig2m = zeros(n,1);

for iter = 1:n

    % beta
    beta = Gen_beta(Y, X, beta_0, b0);  %sig=1
    betam(iter,:) = beta';

    %%sig2
    sig2 = Gen_sig2(Y, X, v_0, d_0, beta);
    sig2m(iter,:) = sig2';

    % Z
     for t = 1:T
        xt = X(t, :)';
        mu = xt'*Beta;

        if Y(t) > 0
            Z(t) = Y(t);
        
        else Y(t) == 0;
            Z(t) = trandn(mu, 1, -100, 0);

        end

    end

end



%Burn-in
betam = betam(n0+1:n,:);
sig2m = sig2m(n0+1:n,1);
MHm = [betam sig2m];

alpha = 0.05;
maxac = 200;
is_plot = 0;
postmom = MHout(MHm,alpha,maxac,is_plot);

% Display
disp('===================================================');
disp('Bayesian estimation');
disp('----------------------------------------------------');
disp(['  Estimate' '     s.e.' '      2.5' '     97.5']);
disp('----------------------------------------------------');
disp([postmom(:,2) postmom(:,3) postmom(:,4) postmom(:,6)]);
disp('----------------------------------------------------');
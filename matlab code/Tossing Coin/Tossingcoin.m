%% Tossing coin
% Sampling Probability of Head
clc;
clear;

%% Step 1 : DGP
% Ym = [1;0;0;0;1;0;0;0;0;0];
T = 1000000;                         % Number of observations (Data information)
theta = 0.9;                           % Probability of Head
Ym = rand(T,1) < theta*ones(T,1);      % Data

%% Step 2 : Frequentist
[T,~] = size(Ym);
NumH = sum(Ym);
Pr_H_Freq = NumH/T;                    % Probability of Head(MLE) = 0.2

%% Step 3 : Bayesian
% Sampling size
n = 10;

% Prior for theta (Prior information)
a0 = 50;
b0 = 50;

% Posterior distribution
% Pr_H_bayes = (NumH + a0)/(a0 + b0 + T);

a1 = a0 + NumH;
b1 = b0 + T - NumH;
Pr_H_Bayes = betarnd(a1,b1,n,1);

%% Step 4 : Display
% Table (This result depends on the prior or data information)
disp('==========================');
disp([' Frequentist','  Bayesian']);
disp('==========================');
% disp ([Pr_H_Freq, Pr_H_Bayes]);
disp([Pr_H_Freq mean(Pr_H_Bayes)]);
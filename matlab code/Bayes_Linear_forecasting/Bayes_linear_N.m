function [betam,sig2m] = Bayes_linear_N(Y,X,beta_0,B_0,v_0,d_0,n0,n1)

% Information
[~,k] = size(X);

% Simulation size
n = n0 + n1;

% pre-allocation
betam = zeros(n,k);
sig2m = zeros(n,1);

% initial value
sig2 = d_0/v_0; % prior mean = d_0/v_0

for iter = 1:n
    
    % Step 1: full conditional posterior dist of beta, given sig2
    %Mean/Variance
    [beta] = Gen_beta(Y,X,B_0,beta_0,sig2);
    betam(iter,:) = beta';

    % Step 2: full conditional posterior dist of sig2, given beta
    [sig2] = Gen_sig2(Y,X,v_0,d_0,beta);
    sig2m(iter,1) = sig2;

end

% Burn-in
betam = betam(n0+1:n,:);
sig2m = sig2m(n0+1:n,1);

end
function [betam,sig2m] = Bayes_linear_N(Y,X,beta_0,B_0,v_0,d_0,n0,n1)

% Information
[T,k] = size(X);

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
    B1 = inv((1/sig2)*X'*X + inv(B_0));
    A = (1/sig2)*X'*Y + inv(B_0)*beta_0;

    % Sampling ans Save
    beta = mvnrnd(B1*A, B1)';
    betam(iter,:) = beta';

    % Step 2: full conditional posterior dist of sig2, given beta
    % Parameters
    v_1 = v_0 + T;
    d_1 = d_0 + (Y-X*beta)'*(Y-X*beta);

    % Sampling and Save
    sig2 = randig(v_1/2,d_1/2,1,1);
    sig2m(iter,1) = sig2;

end

% Burn-in
betam = betam(n0+1:n,:);
sig2m = sig2m(n0+1:n,1);

end
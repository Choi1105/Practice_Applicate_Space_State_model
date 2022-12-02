function [betam,sig2m,yfm,SEm,PPD] = Bayes_linear_N_Pred(Y_Full,X_Full,beta_0,B_0,v_0,d_0,n0,n1,H)

% Information
[~,k] = size(X_Full);
T0 = rows(Y_Full) - H;

% Simulation size
n = n0 + n1;

% pre-allocation for predictive values
PPD = zeros(H,1);
SEm = zeros(H,1);
yfm = zeros(n1,H);

% initial value
sig2 = d_0/v_0; % prior mean = d_0/v_0

for indH = 1:H
    
    % In-sample
    Y = Y_Full(1:T0+indH-1,1);
    X = X_Full(1:T0+indH-1,:);

    % Out-of-sample
    ya = Y_Full(T0+indH,1);
    xf = X_Full(T0+indH,:);

    % Pre-allocation
    betam = zeros(n,k);
    sig2m = zeros(n,1);
    PPDm = zeros(n,1);
    yf = zeros(n,1);

    for iter = 1:n

    % Step 1: full conditional posterior dist of beta, given sig2
    %Mean/Variance
    [beta] = Gen_beta(Y,X,B_0,beta_0,sig2);
    betam(iter,:) = beta';

    % Step 2: full conditional posterior dist of sig2, given beta
    [sig2] = Gen_sig2(Y,X,v_0,d_0,beta);
    sig2m(iter,1) = sig2;

    % Step 3: Forecasting
    yf(iter,1) = mvnrnd(xf*betam(iter,:)',sig2m(iter));

    % Step 4: Compute PPD
    PPDm(iter,1) = normpdf(ya,xf*betam(iter,:)',sqrt(sig2m(iter)));
    end

% Burn-in
betam = betam(n0+1:n,:);
sig2m = sig2m(n0+1:n,1);
PPDm = PPDm(n0+1:n,1);
yf = yf(n0+1:n,1);

% Save
yfm(:,indH) = yf;
SEm(indH) = (ya - meanc(yf))^2;
PPD(indH) = meanc(PPDm);

end

end
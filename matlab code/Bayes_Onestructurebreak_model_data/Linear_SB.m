%% One Structural Break Function

function [beta1m, beta2m, sig21m, sig22m, taum, tau_tsm] = Linear_SB(Y, X, beta_0, B_0, v_0, d_0, a_0, b_0, n0, n1)

% Information
[T, k] = size(X);

% Simulation size
n = n0 + n1;

% Pre-allocation            % 저장할 공간들
taum = zeros(n,1);
tau_tsm = zeros(T,1);
beta1m = zeros(n,k);
beta2m = zeros(n,k);
sig21m = zeros(n,1);
sig22m = zeros(n,1);

% Initial values
sig21 = d_0/v_0;
sig22 = sig21;
tau = round(0.5*T);

for iter = 1:n

    % Step 1: Posterior Conditional Distribution of b, given sig2
    Y1 = Y(1:tau-1);        % t이전 Y1
    X1 = X(1:tau-1,:);      % t이전 X1
    Y2 = Y(tau:end);        % t이후 Y2
    X2 = X(tau:end,:);      % t이후 X2

    beta1 = Gen_beta(Y1,X1,B_0,beta_0,sig21);      
    beta1m(iter,:) = beta1';
    beta2 = Gen_beta(Y2,X2,B_0,beta_0,sig22);
    beta2m(iter,:) = beta2';

    % Step 2: Posterior Conditional Distribution of sig2, given b
    sig21 = Gen_sig2(Y1,X1,v_0,d_0,beta1);
    sig21m(iter,1) = sig21;
    sig22 = Gen_sig2(Y2,X2,v_0,d_0,beta2);
    sig22m(iter,1) = sig22;

    % Step 3: Sampling the breakpoint
    tau = Gen_tau(Y, X, beta1, beta2, sig21, sig22, a_0, b_0);
    taum(iter) = tau;

    if iter > n0
        tau_tsm(tau) = tau_tsm(tau) + 1;
    end

end

% burn-in
beta1m = beta1m(n0+1:n,:);
beta2m = beta2m(n0+1:n,:);
sig21m = sig21m(n0+1:n,:);
sig22m = sig22m(n0+1:n,:);
taum = taum(n0+1:n,:);
tau_tsm = tau_tsm/n1;

end

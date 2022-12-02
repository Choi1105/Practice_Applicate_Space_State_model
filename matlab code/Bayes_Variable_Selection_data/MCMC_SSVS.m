%% SSVS Sampling

function [betam, sig2m, Gam, postmom_beta, postmom_sig2, postmom_Gam] = MCMC_SSVS(Y, X, n0, n1, v_0, d_0, v_00, d_00, v_01, d_01, c0, d0)

%information
n = n0 + n1;
k = cols(X);

% initial values
sig2 = d_0/v_0;
b0 = (d_00/v_00);
b1 = (d_01/v_01);
gam = ones(k, 1);
p = 0.5;

% pre-allocation
Gam = zeros(n, k);
betam = zeros(n, k);
sig2m = zeros(n, 1);

for iter = 1:n

    % Step 1: beta sampling
    beta = Gen_beta_VS(Y, X, sig2, gam, b0, b1);
    betam(iter,:) = beta';

    % Step 2: sig2 sampling
    sig2 = Gen_sig2(Y, X, v_0, d_0, beta);
    sig2m(iter, :) = sig2;

    % Step 3: gam sampling
    gam = Gen_gam(beta, b0, b1, p);
    Gam(iter, :) = gam';

    % Step 4: b0, b1 sampling
    beta_0m = beta(gam == 0);
    b0 = Gen_B(beta_0m, v_00, d_00);

    beta_1m = beta(gam == 1);
    b1 = Gen_B(beta_1m, v_01, d_01);

    % Step 5: p sampling
    p = Gen_p(gam, k, c0, d0);

end

% burn-in
betam = betam(n0+1:n, :);
sig2m = sig2m(n0+1:n, :);
Gam = Gam(n0+1:n, :);

% posterior moments
alpha = 0.05;
maxac = 200;
postmom_beta = MHout(betam, alpha, maxac, 0);
postmom_sig2 = MHout(sig2m, alpha, maxac, 0);
postmom_Gam = MHout(Gam, alpha, maxac, 0);

end
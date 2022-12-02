%% Recursive VAR using Cholesky Decomposition

function [PSIm, MHm] = Recursive_VAR(n0, n1, beta0, B0, nu0, R0, Y, p, mlag)
% Data and Sampling info
k = cols(Y);
n = n0 + n1;

% Make Y and X
[Ym, YLm] = makeYX(Y,p);

% Pre-allocation for beta, Omega, IR
betam = zeros(n1,p*k*k);
Omegam = zeros(n1, k^2);
PSIm = zeros(n1, mlag+1, k^2);

% Initial values
Phi = reshape(beta0, p*k, k);
Omega_inv = nu0*R0;

% MCMC START!!
for iter = 1:n

    % Phi Sampling
    [Phi, Fm, beta] = Gen_Phi(Ym, YLm, Phi, beta0, B0, Omega_inv);

    % Omega Sampling
    [Omega, Omega_inv] = Gen_Omega(Ym, YLm, beta, nu0, R0);

    % After burn-in, Save
    if iter > n0
        PSIm = Gen_ImRes(Omega, Fm, mlag, n0, PSIm, iter);
        betam(iter-n0,:) = beta';
        Omegam(iter-n0,:) = vec(Omega)';
    end

end

% Joint poesterior samples
MHm = [betam Omegam];

end
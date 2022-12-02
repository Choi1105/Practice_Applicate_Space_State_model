%% Sampling Omega

function [Omega, Omega_inv] = Gen_Omega (Y, X, beta, nu0, R0)

% Sample size(T), # of variables(k)
T = rows(Y);
k = cols(Y);

% Sampling full conditional distribution
% Step 1: Compute ehat^2
ehat2 = zeros(k,k);
for t = 1:T
    Xt = X(:,:,t);
    ehat = Y(t,:)' - Xt*beta;
    ehat2 = ehat2 + ehat*ehat';
end

% Step 2: Compute Inverse wishart parameters
R1 = R0 + ehat2;
nu1 = nu0 + T;

% Step 3: Sampling Omega
Omega = iwishrnd(R1, nu1);
Omega_inv = inv(Omega);

end
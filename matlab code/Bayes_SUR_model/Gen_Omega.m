%% Sampling Omega
% Input: Data(Y,X), Prior(nu0,R0), Given values(beta)
% Output: Sampled values(Omega,Omega_inv)

function [Omega,Omega_inv] = Gen_Omega(Y,X,beta,nu0,R0)

% Sample size(T)
T = rows(Y);

% Sampling full conditional distribution
% STEP 1: Compute ehat^2
ehat2 = 0;
for t = 1:T
    Xt = X(:,:,t);
    ehat = Y(t,:)' - Xt*beta;
    ehat2 = ehat2 + ehat*ehat';
end

% STEP 2: Compute Inverse wishart parameters
R1_inv = ehat2 + invpd(R0);
R1 = invpd(R1_inv);
nu1 = T + nu0;

% STEP 3: Sampling Omega
Omega_inv = randwishart(R1,nu1);
Omega = invpd(Omega_inv);

end
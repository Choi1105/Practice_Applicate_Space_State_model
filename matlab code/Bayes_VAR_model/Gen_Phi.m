%% Sampling Phi

function [Phi, Fm, beta] = Gen_Phi (Y, X, Phi0, beta0, B0, Omega_inv)

% Sample size(T), # of variables(k), lag(p)
[k, pkk, T] = size(X);
p = pkk/k^2;

% Samlping Full conditional distribution
% Step 1: Compute XX and XY
XX = 0;
XY = 0;
for t = 1:T
    Xt = X(:,:,t);
    XX = XX + Xt'*Omega_inv*Xt;
    XY = XY + Xt'*Omega_inv*Y(t,:)';
end

% Step 2: Compute full conditional mean/variance
B1 = inv(XX + inv(B0));
A = XY + inv(B0)*beta0;

% Step 3: beta sampling
beta = mvnrnd(B1*A, B1)';

% Make F matrix, p*k by p*k
Phi = reshape(beta, p*k, k);                    % p*k by k
Fm = [Phi'; eye((p-1)*k), zeros(k*(p-1), k)];   % p*k by p*k

% Step 4: Check stationary condition
% If eigenvalues > 1, Reject (Retain previous values)
% Otherwise         , Accept (Save new values)
eigF = eig(Fm);
if maxc(abs(eigF)) >= 1
    Phi = Phi0;
    Fm = [Phi'; eye((p-1)*k), zeros(k*(p-1), k)];
    beta = vec(Phi);
end

end
%% Compute Unconditional Variance of Factors
% Assume Var[f(t)] = Var[f(t-1)]
% Var[f(t)] = reshape(NewQ,k,k)
% WHERE NewQ = INV[eye(k^2) - kron(G)]*vec(Q)

% Input: F, Omega (Model parameters)
% Output: P_ll (Unconditional Variance)

function P_ll = makeP_ll(F,Omega) 

k = rows(F); 
k2 = k^2;
G2 = kron(F,F); 
eyeG2 = eye(k2) - G2; 
omegavec = reshape(Omega,k2,1);
P_ll = (eyeG2)\omegavec; 
P_ll = reshape(P_ll,k,k)';
P_ll = (P_ll + P_ll')/2;

end    
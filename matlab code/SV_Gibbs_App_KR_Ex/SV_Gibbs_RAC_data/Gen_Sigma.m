%% Generate Sig2
% input: hm(Data), v_,d_(prior), mu_phi(other parameters)
% output: sig2, sig2_inv(sampled values)

function [sig2_inv,sig2] = Gen_Sigma(Y,X,v0,d0,mu_phi)

% Sample size
T = rows(Y);

% residual term
ehat = Y - X*mu_phi; 

% Full conditional distribution
v1 = v0 + T;
d1 = d0 + ehat'*ehat;

% Sampling from IG
sig2 = randig(v1/2,d1/2,1,1);
sig2_inv = 1/sig2;

end
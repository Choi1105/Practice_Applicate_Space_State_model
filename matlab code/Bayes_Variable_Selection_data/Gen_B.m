%% B Sampling

function B = Gen_B(beta, v0, d0)

% parameters
d1 = d0 + beta'*beta;
v1 = v0 + rows(beta);

% Sampling
B = randig(v1/2, d1/2, 1, 1);

end
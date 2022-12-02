%% p Sampling

function p = Gen_p(gam, k, c0, d0)

% parameters
k1 = sumc(gam);
k0 = k - k1;

% sampling
p = betarnd(c0 + k1, d0 + k0);

end
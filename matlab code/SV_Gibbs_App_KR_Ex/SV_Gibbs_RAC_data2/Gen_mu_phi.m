%% Generate Mu and Phi
% input: hm(Data), b_,precb_(prior),
%        sig2_inv(sig2 inverse),mu_phi0(previous value)
% output: mu_phi (sampled value)

function [mu_phi] = Gen_mu_phi(Y,X,b_,precb_,sig2_inv,mu_phi0)

%  # of parameters
k = cols(X);                          

% Compute full conditional variance/mean(FCV,FCM)
XX = X'*X;
XY = X'*Y;
varb1 = invpd(precb_ + sig2_inv*XX);  % FCV
b1 = varb1*(precb_*b_ + sig2_inv*XY); % FCM

% Beta sampling
mu_phi = b1 + chol(varb1)'*randn(k,1);

% Stationary conditions
% Hold    : Save sampled value
% NOT Hold: Save previous value
if abs(mu_phi(2)) > 1
    mu_phi = mu_phi0;
end

end
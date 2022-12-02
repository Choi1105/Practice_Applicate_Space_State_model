%% Sampling Beta
% Input: Data(Y,X), Prior(b0,var0), Given values(Phi0,Omega_inv)
% Output: Sampled values(Phi,Fm,beta)

function beta = Gen_Beta(Y,X,b0,precb0,Omeg_inv)

% Sample size(T), # of total coefficients(k)
T = rows(Y);
k = rows(b0);

% Sampling full conditional distribution
% STEP 1: Compute XX and XY
XX = 0;
XY = 0;
for t = 1:T
    Xt = X(:,:,t);
    XX = XX + Xt'*Omeg_inv*Xt;
    XY = XY + Xt'*Omeg_inv*Y(t,:)';
end

% STEP 2: compute full conditional mean/variance
% mean
B1_inv = precb0 + XX;
B1_inv = 0.5*(B1_inv + B1_inv');
B1 = invpd(B1_inv);
B1 = 0.5*(B1 + B1');
A = XY + precb0*b0; 
BA = B1*A; 
% variance
Chol_B1 = cholmod(B1)';

% STEP 3: beta sampling
beta = BA + Chol_B1*randn(k,1); 

end
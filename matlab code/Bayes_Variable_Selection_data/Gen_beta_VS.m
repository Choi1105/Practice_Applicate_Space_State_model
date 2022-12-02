%% Beta Sampling

function [beta] = Gen_beta_VS(Y, X, sig2, gam, b0, b1)

% Prior
B_0 = diag((1 - gam)*b0 + gam*b1);

% Mean/Variance
b1 = inv((1/sig2)*X'*X + inv(B_0));
A = (1/sig2)*X'*Y;

    % Sampling ans Save
    beta = mvnrnd(b1*A, b1)';
    
end
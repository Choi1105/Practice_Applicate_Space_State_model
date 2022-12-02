%% Beta Sampling Function

function beta = Gen_beta(Y,X,beta_0, B_0,sig2)
    
% Mean/Variance
B1 = inv((1/sig2)*X'*X + inv(B_0));
A = (1/sig2)*X'*Y + inv(B_0)*beta_0;

% Sampling ans Save
beta = mvnrnd(B1*A, B1)';
end
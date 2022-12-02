function beta = Gen_beta(Y, X, beta_0, B_0)

% Mean/Variance
B1 = inv(X'*X + inv(B_0));       % 2 by 2
A = X'*Y + inv(B_0)*beta_0;      % 2 by 1

% Sampling
beta = mvnrnd(B1*A,B1)';

end
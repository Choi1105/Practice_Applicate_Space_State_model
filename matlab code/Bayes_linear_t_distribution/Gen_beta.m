function [beta] = Gen_beta(Y,X,B_0,beta_0,sig2,lamda)
    lamda = diag(1);
    B1 = inv((1/sig2)*X'*lamda*X + inv(B_0));
    A = (1/sig2)*X'*lamda*Y + inv(B_0)*beta_0;

    % Sampling ans Save
    beta = mvnrnd(B1*A, B1)';
end
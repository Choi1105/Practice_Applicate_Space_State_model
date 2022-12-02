function [sig2] = Gen_sig2(Y,X,a_0,d_0,beta,lamda)
% Parameters
    lamda = diag(1);
    a_1 = a_0 + rows(Y);
    d_1 = d_0 + (Y-X*beta)'*lamda*(Y-X*beta);

    % Sampling and Save
    sig2 = randig(a_1/2,d_1/2,1,1);
end
function log_density=MVN_PDF_LOG(X,mu,sigma)

% Purpose: 
% Log density (p.d.f) of multivariate normal distribution
% -----------------------------------
% Density:
% f(x) = (2*pi)^(-k/2) * det(cov_mat)^(-1/2) * exp(-(x-mu)'*inv(cov_mat)*(x-mu)/2)
% where k is the dimension of MVN
% E(X) = mu, Cov(X) = cov_mat
% -----------------------------------
% Algorithm: 
% Cholesky decomposition of covariance matrix. ¦² = P*P¡¯
% det(¦²) = prod(diag(P))
% (X-¦Ì)' * inv(¦²) * (X- ¦Ì) = [inv(P) * (X-¦Ì)]¡¯ *  [inv(P) * (X-¦Ì)]
% -----------------------------------
% Usage:
% X = points of evaluation (k*1)
% mu = mean vector of multivariate normal distribution (k*1)
% sigma = covariance matrix of multivariate normal distribution (k*k)
% -----------------------------------
% Returns:
% log_density = log density evaluated at points X
% -----------------------------------
% Notes:
% Support matrix input of X, mu, which should be k * n
% It will return a vector density with conformable size.
%
% Written by Hang Qian, Iowa State University
% Contact me:  matlabist@gmail.com


% smart adjustment of matrix dimension
tranpose_flag = 0;
dim = size(sigma,1);
[nrow1,ncol1] = size(X);
if (nrow1 ~= dim) &&  (ncol1 == dim)
    X = X';
    [nrow1,ncol1] = size(X);
    tranpose_flag = 1;
end

[nrow2,ncol2] = size(mu);
if (nrow2 ~= dim) &&  (ncol2 == dim)
    mu = mu';
    [nrow2,ncol2] = size(mu);
    if nrow1 == 1
        tranpose_flag = 1;
    end
end

nobs = max(ncol1,ncol2);
if (nobs > ncol1) && (nrow1 > 1) 
    X = X(:,ones(nobs,1));
end

if (nobs > ncol2) && (nrow2 > 1) 
    mu = mu(:,ones(nobs,1));
end

% Compute density
constant = -0.5*dim*log(2*pi);

P = chol(sigma,'lower');
det_term = -sum(log(diag(P)));

X_trans = P\(X-mu);
exp_term = -0.5 * sum(X_trans.*X_trans);

log_density = constant + det_term + exp_term;

if tranpose_flag == 1
    log_density = log_density';
end
    
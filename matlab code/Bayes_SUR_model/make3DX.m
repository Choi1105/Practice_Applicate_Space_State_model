%% Make X(3-dim)
% Input: Data(Y), number of regression(p)
% Output: 3-dim Exogenous variables

function Xm = make3DX(X,p)

% Data information
T = rows(X);    % sample size
k = cols(X);    % total regression coefficient
ki = k/p;       % # of regression coefficient in each regression

% 3-dimension, p by k by T
Xm = zeros(p,k,T);
for t = 1:T
    xt = zeros(p,k);
    for indp = 1:p
        xt(indp,(indp-1)*ki+1:indp*ki) = X(t,(indp-1)*ki+1:indp*ki);
    end
    Xm(:,:,t) = xt; 
end

end
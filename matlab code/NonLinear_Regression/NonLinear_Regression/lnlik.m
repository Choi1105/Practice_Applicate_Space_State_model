%% Compute log likelihood function

function [lnL] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:end);

% Parameter
rho = theta(1); 
ibar = theta(2);
alpha = theta(3);
beta = theta(4); 
sig2 = theta(5);

% Compute log density for t=1~T
lnLm = zeros(T,1);
for t = 1:T
    yt = Y(t);
    xt = X(t,:);
    yf = (1-rho)*ibar*xt(1) + (1-rho)*alpha*xt(2) + (1-rho)*beta*xt(3) + rho*xt(4);
    lnLm(t) = lnpdfn(yt, yf, sig2);
end

% Sum
lnL = sumc(lnLm);

end
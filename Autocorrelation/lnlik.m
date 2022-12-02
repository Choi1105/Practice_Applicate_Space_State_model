%% Compute log likelihood function

function [lnL] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:end);

% Parameter
beta1 = theta(1);
beta2 = theta(2); 
rho = theta(3);
sig2 = theta(4); 

% Compute log density for t=1~T
lnLm = zeros(T,1);
for t = 2:T
    yt = Y(t);
    yf = beta1 + beta2*X(t,2) ...
       + rho*(Y(t-1) - beta1 - beta2*X(t-1,2));
    % yf = (1 - rho)*beta1*X(t,1) + rho*Y(t-1) + beta2*(X(t,2) - rho*X(t-1,2));
    lnLm(t) = lnpdfn(yt, yf, sig2);
end

% Sum
lnL = sumc(lnLm);

end
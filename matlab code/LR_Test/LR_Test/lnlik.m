%% Compute log likelihood function

function [lnL] = lnlik(theta,Data)

% Data
T = rows(Data);

Y = Data(:,1);
X = Data(:,2:end);

% Parameter
beta = theta(1:2);
sig2 = theta(3);

% Compute log density for t=1~T
lnLm = zeros(T,1);
for t = 1:T
    yt = Y(t);
    xt = X(t,:);
    yf = xt*beta;
    lnLm(t) = lnpdfn(yt, yf, sig2);
end

% Sum
lnL = sumc(lnLm);

end
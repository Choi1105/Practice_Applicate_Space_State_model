%% Compute log likelihood function

function [lnL] = lnlik(theta,Data)

% Data
T = rows(Data);
Y = Data(:,1);
X = Data(:,2:end-1);
Z = Data(:,end);

% Parameter
beta = theta(1:2); 
a0 = theta(3);
a1 = theta(4); 

% Compute log density for t=1~T
lnLm = zeros(T,1);
for t = 1:T
    yt = Y(t);
    xt = X(t,:);
    zt = Z(t);
    yf = xt*beta;
    sig2t = a0 + a1*zt;
    lnLm(t) = lnpdfn(yt, yf, sig2t);
end

% Sum
lnL = sumc(lnLm);

end
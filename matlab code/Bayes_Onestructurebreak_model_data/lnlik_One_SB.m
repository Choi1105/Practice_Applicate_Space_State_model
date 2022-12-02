%% Compute lnlikelihood

function [lnL] = lnlik_One_SB(X, Y, SB, beta1, sig21, beta2, sig22)

T = rows(Y);
lnLm = zeros(T,1);

for t = 1:T
    
    if t < SB
        lnf = log(normpdf(Y(t), X(t,:)*beta1, sqrt(sig21)));
    else
        lnf = log(normpdf(Y(t), X(t,:)*beta2, sqrt(sig22)));
    end
    lnLm(t) = lnf;
end

lnL = sumc(lnLm);

end
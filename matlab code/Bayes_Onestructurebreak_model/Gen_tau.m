%% Tau Sampling Function

function [tau] = Gen_tau(Y, X, beta1, beta2, sig21, sig22, a_0, b_0)

T = rows(X);
Probm = zeros(T,1);

for t = a_0+1:b_0-1

    SB = t;
    Probm(t) = exp(lnlik_One_SB(X, Y, SB, beta1, sig21, beta2, sig22));

end

tau_prb = Probm/(sumc(Probm));
tau = discret1 (tau_prb,1);

end
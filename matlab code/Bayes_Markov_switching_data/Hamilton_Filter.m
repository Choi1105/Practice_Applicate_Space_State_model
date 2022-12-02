%% Hamilton Filter: Forward Recursion
% To compute Filtered probability

function Filtpm = Hamilton_Filter(Y, X, beta1, beta2, sig21, sig22, p11, p22)
% Data information
T = rows(Y);  % Sample size
ns = 2;       % Regimes

% Unconditional Prob.
prb1_LL = (1 - p11) / (2 - p11 - p22);
prb2_LL = 1 - prb1_LL;

% Pre-allocation for Filtered Prob.
Filtpm = zeros(T, ns);

for t = 1:T

    % Pr[S(t),S(t-1)|I(t-1)]
    prb11_tL = prb1_LL*p11;
    prb12_tL = prb1_LL*(1-p11);
    prb21_tL = prb2_LL*(1-p22);
    prb22_tL = prb2_LL*p22;

    % Pr[S(t)|I(t-1)]
    prb1_tL = prb11_tL + prb21_tL;
    prb2_tL = prb12_tL + prb22_tL;

    % f[y(t)|I(t-1)]
    pdf1 = mvnpdf(Y(t), X(t,:)*beta1, sig21);
    pdf2 = mvnpdf(Y(t), X(t,:)*beta2, sig22);

    pdf = pdf1*prb1_tL + pdf2*prb2_tL;

    % Pr[S(t)|I(t)]
    prb1_tt = (pdf1*prb1_tL)/pdf;
    prb2_tt = (pdf2*prb2_tL)/pdf;

    % Save Filtered Prob
    Filtpm(t,:) = [prb1_tt, prb2_tt];

    % Updating
    prb1_LL = prb1_tt;
    prb2_LL = prb2_tt;

end
%% Generate P from Beta distribution

function [p11, p22, P] = Gen_Prb(S, a10, b10, a20, b20)

% # of regimes
ns = maxc(S);

% Count # of regime switches
Nm = CountTransitions(S);

% Sampling p11, p22
a11 = Nm(1,1) + a10;
b11 = Nm(1,2) + b10;
p11 = betarnd(a11, b11);

b21 = Nm(2,1) + b20;
a21 = Nm(2,2) + a20;
p22 = betarnd(a21, b21);

% Transition Probability Matrix
P = zeros(ns, ns);
P(1,1) = p11;
P(1,2) = 1 - p11;
P(2,2) = p22;
P(2,1) = 1 - p22;

end
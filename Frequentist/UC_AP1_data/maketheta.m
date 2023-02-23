%% Transform parameters
% input: psi
% output: theta (transfomred parameters)
% This code is necessary when we need some parameter constraints

function [theta] = maketheta(psi, Sn)

% Variance in Transtion equation > 0
theta = psi;
theta(Sn.indQ) = exp(psi(Sn.indQ));

end
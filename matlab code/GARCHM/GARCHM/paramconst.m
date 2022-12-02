%% Parameter constraint
% input: psi, Sn (paramters, structure parameters)
% output: valid = 1, parameter constraint satisfied
%         valid = 0, otherwise

function [valid] = paramconst(theta,Data)
% Pre-allocation
validm = ones(30,1);

% a0 > 0
validm(1) = theta(3) > 0;

% a1 + gam1 < 1
validm(2) = theta(4) + theta(5) < 1;

% a1 > 0
validm(3) = theta(4) > 0;

% gam > 0
validm(4) = theta(5) > 0;

valid = minc(validm);

end
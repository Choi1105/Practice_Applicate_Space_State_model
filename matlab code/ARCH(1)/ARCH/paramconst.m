%% Parameter constraint
% input: psi, Sn (paramters, structure parameters)
% output: valid = 1, parameter constraint satisfied
%         valid = 0, otherwise

function [valid] = paramconst(theta,Data)

% Pre-allocation
validm = ones(30,1);

% a0 > 0
validm(1) = theta(2) > 0;

% 0 < a1 < 1
validm(2) = theta(3) > 0;
validm(3) = theta(3) < 1;

valid = minc(validm);

end
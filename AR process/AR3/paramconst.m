%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% sig2 > 0 
validm(1) = theta(5) > 0;

% stationary
validm(2) = abs(theta(2)+theta(3)+theta(4)) < 1;

valid = min(validm);

end


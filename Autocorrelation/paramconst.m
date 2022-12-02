%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% -1 < rho < 1
validm(1) = theta(3) > -1; 
validm(2) = theta(3) < 1; 

% sig2 > 0
validm(3) = theta(4) > 0; 

valid = min(validm); 

end

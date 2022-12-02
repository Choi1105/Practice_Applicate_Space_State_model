%% parameter constraints

function [valid] = paramconst(theta,Data)

validm = ones(10,1);

% Sig2 > 0
validm(1) = theta(5) > 0; 

% 0 < rho < 1
validm(2) = theta(1) > 0; 
validm(3) = theta(1) < 1; 

valid = min(validm); 

end

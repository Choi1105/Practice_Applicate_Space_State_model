%% Parameter constraint
% input: psi, Sn (paramters, structure parameters)
% output: valid = 1, parameter constraint satisfied
%         valid = 0, otherwise

function [valid] = paramconst(psi,Sn)
% Define structural parameters
indQ = Sn.indQ;
indR = Sn.indR;

% Pre-allocation
validm = ones(30,1);

% Except for "NAN" parameters 
if maxc(isnan(psi)') == 1
    validm(1) = 0;
end 

% Variance > 0
validm(2) = minc(psi(indR)) > -8;
validm(3) = minc(psi(indQ)) > -8;

% Check whether all constranit satisfied
valid = minc(validm); 


end
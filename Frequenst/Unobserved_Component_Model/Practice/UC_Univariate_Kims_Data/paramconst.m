%% Parameter constraint
% input: psi, Sn (paramters, structure parameters)
% output: valid = 1, parameter constraint satisfied
%         valid = 0, otherwise

function [valid] = paramconst(psi,Sn)

% Define structural parameters
indQ = Sn.indQ;
indF = Sn.indF;

% Pre-allocation
validm = ones(30,1);

% Except for "NAN" parameters 
if maxc(isnan(psi)') == 1
    validm(1) = 0;
end 

% Transition equation variance > 0  % control variance
%validm(2) = minc(psi(indQ)) > -5; % control variance

% If sig2v and sig2w and sig2e is too low, then Function like Z(t) = Z(t-1)
% so you have to constrain the parameter, especially variance of Error term

% Stationarity condition
Phim = psi(indF);
F = [Phim(1) Phim(2);1 0];
eigF = eig(F);
validm(3) = maxc(abs(eigF)) < 1;

% Check whether all constranit satisfied
valid = minc(validm); 

end
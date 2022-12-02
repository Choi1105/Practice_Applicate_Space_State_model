%% Compute Unconditional Mean of Factors
% Assume E[f(t)] = E[f(t-1)]
% E[f(t)] = INV(I(k) - F)*mu

% Input: F, Mu (Model parameters)
% Output: f_ll (Unconditional Mean)

function f_ll = makef_ll(F,Mu)
        
k = rows(F);    
eyeF = eye(k)-F; 
f_ll = eyeF\Mu;

end
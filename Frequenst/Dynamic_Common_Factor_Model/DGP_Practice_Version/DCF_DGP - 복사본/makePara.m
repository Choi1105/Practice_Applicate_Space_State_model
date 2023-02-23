%% Make State Space Parameter Form
% [Input]
% theta: Set of parameters
% Sn: Structure parameters

% [Output]
% Mu,F,Q,H,C,R: Set of State Space Parameter

function [C,H,R,Mu,F,Q] = makePara(theta, Sn) 

% Structure parameters
indH = Sn.indH;
indR = Sn.indR;
indMu = Sn.indMu;
indF = Sn.indF;
indQ = Sn.indQ;

% Transition equation
Mu = theta(indMu);

F = theta(indF);

Q = theta(indQ);

% Measurement equation
C = zeros(2,1);

H = ones(2,1); 
H(2) = theta(indH); 

R = diag(theta(indR)); 

end 
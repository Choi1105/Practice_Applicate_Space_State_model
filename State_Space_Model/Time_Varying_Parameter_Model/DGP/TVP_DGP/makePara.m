%% Make State Space Parameter Form
% [Input]
% theta: Set of parameters
% Sn: Structure parameters

% [Output]
% Mu,F,Q,H,C,R: Set of State Space Parameter

function [C,H,R,Mu,F,Q] = makePara(theta, Sn) 

% Structure parameters
xm = Sn.xm;
indR = Sn.indR;
indQ = Sn.indQ;

% Transition equation
Mu = 0;

F = 1;

Q = theta(indQ); 

% Measurement equation
C = 0;

H = xm;

R = theta(indR);

end 
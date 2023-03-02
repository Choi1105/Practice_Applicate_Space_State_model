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
Mu = [0 0 0 0]';

F = eye(4);

sig2v = theta(indQ);
Q = diag([sig2v(1),sig2v(2),sig2v(3),sig2v(4)]);

% Measurement equation
C = 0;

H = xm;

R = theta(indR);

end 
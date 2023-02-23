%% Make State Space Parameter Form
% [Input]
% theta: Set of parameters
% Sn: Structure parameters

% [Output]
% Mu,F,Q,H,C,R: Set of State Space Parameter

function [C,H,R,Mu,F,Q] = makePara(theta, Sn) 

% Structure parameters
indMu = Sn.indMu;
indF = Sn.indF;
indQ = Sn.indQ;

% Transition equation
Mu = [theta(indMu);0;0];

Phim = theta(indF);
F = [ 1 0 0 ; 0 Phim(1) Phim(2) ; 0 1 0 ];

Q = diag([theta(indQ);0]); 

% Measurement equation
C = 0;

H = [1 1 0];

R = 0;

end 
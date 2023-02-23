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
indH = Sn.indH;
indR = Sn.indR;
indQ = Sn.indQ;

% Transition equation
Mu = [0;0;0;0;0;0];

Phim = theta(indF);
F = [1 1 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0 ; 0 0 0 Phim(1) Phim(2) 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0];

Q = diag([theta(indQ);0;0]); 


% Measurement equation
C = 0;
Am = indH;
H = [1 0 0 1 0 0 ; 0 0 1 Am(1) Am(2) Am(3)];

R = diag([0;theta(indR)]);

end 
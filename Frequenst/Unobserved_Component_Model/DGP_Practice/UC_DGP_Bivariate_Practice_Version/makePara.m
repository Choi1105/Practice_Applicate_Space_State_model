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
indR = Sn.indR;
indH = Sn.indH;

% Transition equation
Mu = [0;0;0;0;0;0];

Phim = theta(indF);
F = [1 1 0 0 0 0 ; 0 1 0 0 0 0 ;  0 0 1 0 0 0  ; 0 0 0 Phim(1) Phim(2) 0 ; 0 0 0 1 0 0 ; 0 0 0 0 1 0];

Q = diag([theta(indQ);0;0]); 

% Measurement equation
C = 0;

alp = indH;
H = [1 0 0 1 0 0 ; 0 0 1 alp(1) alp(2) alp(3)];

R = [0 0; theta(indR) 0 ];

end 
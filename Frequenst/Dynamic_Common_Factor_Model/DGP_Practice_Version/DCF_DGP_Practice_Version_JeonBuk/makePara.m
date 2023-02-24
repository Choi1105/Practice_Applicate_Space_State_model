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
Mu = [theta(indMu) 0 0 ]';

phim = theta(indF);
F = [phim(1) phim(2) 0 ; 1 0 0 ; 0 1 0];

Q = diag([theta(indQ);0;0]);

% Measurement equation
C = zeros(7,1);

Gam = theta(indH); 
H = [1 0 0 ; Gam(1) 0 0 ; Gam(2) 0 0 ; Gam(3) 0 0 ; Gam(4) 0 0 ; Gam(5) 0 0 ; Gam(6) 0 0 ];


R = diag(theta(indR));

end 
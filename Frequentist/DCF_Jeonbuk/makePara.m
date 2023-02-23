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
Mu = [theta(indMu);0];

phim = theta(indF);
F = [ phim(1) phim(2) ; 1 0 ];

Q = diag([theta(indQ);0]);

% Measurement equation
C = zeros(7,1);

gam = theta(indH); 
H = [1 0 ; gam(1) 0 ; gam(2) 0 ; gam(3) 0 ; gam(4) 0 ; gam(5) 0 ; gam(6) 0 ];

sig = theta(indR);
R = diag([sig(1);sig(2);sig(3);sig(4);sig(5);sig(6);sig(7)]); 

end 
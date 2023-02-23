%% Make State Space Parameter Form
% [Input]
% theta: Set of parameters
% Sn: Structure parameters

% [Output]
% Mu,F,Q,H,C,R: Set of State Space Parameter

function [C,H,R,Mu,F,Q] = makePara(theta, Sn) 

% Structure parameters
indR = Sn.indR;
indMu = Sn.indMu;
indF = Sn.indF;
indQ = Sn.indQ;

% Transition equation
MMu = theta(indMu);
Mu = [MMu(1); MMu(2) ; MMu(3)];

phim = theta(indF);
F = [phim(1) 0 0 ; 0 phim(2) 0 ; 0 0 phim(3)];

sigv = theta(indQ);
Q = [sigv(1) sigv(4) sigv(5) ; sigv(4) sigv(2) sigv(6) ; sigv(5) sigv(6) sigv(3)];

% Measurement equation
C = 0;

H = readmatrix("Lambda.xlsx");

R = theta(indR); 

end 
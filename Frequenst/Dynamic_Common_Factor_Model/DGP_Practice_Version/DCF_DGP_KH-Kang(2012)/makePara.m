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
M = theta(indMu);
Mu = [M(1) ; M(2); M(3)];

phi = theta(indF);
F = diag([phi(1);phi(2);phi(3)]);

sigv = theta(indQ);
Q = [sigv(1) sigv(4) sigv(5) ; sigv(4) sigv(2) sigv(6) ; sigv(5) sigv(6) sigv(3)];

% Measurement equation
C = zeros(10,1);

tau = [3 6 9 12 18 24 30 36 60 120];
lambda = 0.15;
Lambda_M = makeLambda(lambda,tau);
H = Lambda_M;

R = diag(theta(indR)); 

end 
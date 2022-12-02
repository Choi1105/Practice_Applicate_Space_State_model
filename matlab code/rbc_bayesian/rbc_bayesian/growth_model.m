% GROWTH_MODEL.M
% This program describes the model 

function [f,x,xp,y,yp,logTransY,logTransX] = growth_model(unitFree,logApprox)


%Define parameters
% syms BETTA GAMA DELTA ALFA PSSI RHOA RHOG OMEGA THETA PHI RSTAR DBAR G 
syms BETTA GAMA DELTA ALFA PSSI RHOA RHOG OMEGA THETA PHI RSTAR DY G 

%Define variables 
syms c cp k kp k1 k1p a ap h hp d dp  yy yyp invs invsp tb tbp la lap tby tbyp cay cayp gy gcc giv yyback cback ivback gyp gccp givp yybackp cbackp ivbackp g gp gback gbackp

%Give functional form for utility, production, and interest-rate premium functions

%argument in expressions related to capital adjustment costs
ac = kp/k*g-G;
acp = k1p/k1*gp-G;

%Trade balance
tb = yy - c - invs - PHI/2 * (kp/k*g -G)^2*k;

%bp = yyp - cp - invsp - PHI/2 * ac^2;



%Interest Rate
% r = RSTAR + PSSI * (exp(d-DBAR) - 1);
% syms YY
r = RSTAR + PSSI * (exp(dp/yy-DY) - 1); 

%Martinal utility of consumption
la = (c - THETA/OMEGA * h^OMEGA)^(-GAMA);

lap = (cp - THETA/OMEGA * hp^OMEGA)^(-GAMA);

% la = (c - 1/OMEGA * h^OMEGA)^(-GAMA);
% 
% lap = (cp - 1/OMEGA * hp^OMEGA)^(-GAMA);

%Write equations (e1, e2,...en)
%Note: we take a linear, rather than log-linear, approximation with respect to tb, the trade balance)
e1 = d - tb - dp*g/r;

e2 = -yy + a * k^ALFA *  (g*h)^(1-ALFA);

e3 = -invs + kp*g - (1-DELTA) *k;

e4 = - la  + BETTA / g^GAMA * r * lap;

e5 = -THETA* h^(OMEGA-1) + (1-ALFA) * a * g^(1-ALFA) * (k/h)^ALFA;
% e5 = - h^(OMEGA-1) + (1-ALFA) * a * g^(1-ALFA) * (k/h)^ALFA;

e6 = -la * (1+ PHI * ac) + BETTA/ g^GAMA * lap * (1 - DELTA + ALFA * ap * (gp*hp/kp)^(1-ALFA) + PHI * k1p/k1*gp * acp - PHI/2 * acp^2 );

e7 = -k1 + kp;

e8 = -log(ap) + RHOA * log(a); 

e9 = -log(gp/G) + RHOG * log(g/G); 

e10 = -tby + tb / yy; 

e11 = -gy + yy/yyback*gback;

e12 = -gcc + c/cback*gback;

e13 = -giv + invs/ivback*gback;

e14 = -yybackp + yy;

e15 = -cbackp + c;

e16 = -ivbackp + invs;

e17 = -gbackp + g;

%Create function f
f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17];

% Define the vector of controls, y, and states, x
x = [yyback cback ivback gback k d g a];
xp = [yybackp cbackp ivbackp gbackp kp dp gp ap];

y = [gy gcc giv tby h yy c invs k1];
yp = [gyp gccp givp tbyp hp yyp cp invsp k1p];

%variables to substitute from levels to logs
logvar = [gy gcc giv h yy c invs k1 yyback cback ivback gback k d g a];
logvarp = [gyp gccp givp hp yyp cp invsp k1p yybackp cbackp ivbackp gbackp kp dp gp ap];

lv = [logvar logvarp];

% For the log-approximation: Make f a function of the logarithm of the state 
% and control vector
if logApprox == 1
    f = subs(f, lv, exp(lv));
    logTransX = ones(1,8);
    logTransY = ones(1,8);
else
    % The indicator variables
    logTransX = zeros(1,8);
    logTransY = zeros(1,8);
end
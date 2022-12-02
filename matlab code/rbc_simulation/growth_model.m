% GROWTH_MODEL.M
% This program describes the model 

function [f,x,xp,y,yp,logTransY,logTransX] = growth_model(unitFree,logApprox)

%Define parameters
% syms BETTA GAMA DELTA ALFA PSSI RHOA RHOG OMEGA THETA PHI RSTAR DBAR G 
syms BETTA GAMA DELTA ALFA PSSI RHOA RHOG OMEGA THETA RSTAR DY YY

%Define variables 
syms c cp k kp k1 k1p a ap h hp d dp yy yyp invs invsp tb tbp la lap tby tbyp cay cayp gy gcc giv yyback cback ivback gyp gccp givp yybackp cbackp ivbackp

%Give functional form for utility, production, and interest-rate premium functions

%argument in expressions related to capital adjustment costs

%Trade balance
tb = yy - c - invs;

%Interest Rate
% r = RSTAR + PSSI * (exp(d-DBAR) - 1);
% syms YY
r = RSTAR + PSSI * (exp(dp/YY-DY) - 1); 

%Martinal utility of consumption
la = (c - THETA/OMEGA * h^OMEGA)^(-GAMA);
lap = (cp - THETA/OMEGA * hp^OMEGA)^(-GAMA);

%Write equations (e1, e2,...en)
%Note: we take a linear, rather than log-linear, approximation with respect to tb, the trade balance)
e1 = d - tb - dp/r;

e2 = -yy + a * k^ALFA * h^(1-ALFA);

e3 = -invs + kp - (1-DELTA) *k;

e4 = - la  + BETTA * r * lap;

e5 = -THETA* h^(OMEGA-1) + (1-ALFA) * a * (k/h)^ALFA;

e6 = -la + BETTA * lap * (1 - DELTA + ALFA * ap * (hp/kp)^(1-ALFA));

e7 = -k1 + kp;

e8 = -log(ap) + RHOA * log(a); 

e9 = -tby + tb / yy; 

e10 = -gy + yy/yyback;

e11 = -gcc + c/cback;

e12 = -giv + invs/ivback;

e13 = -yybackp + yy;

e14 = -cbackp + c;

e15 = -ivbackp + invs;

%Create function f
f = [e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15];

% Define the vector of controls, y, and states, x
x = [yyback cback ivback k d a];
xp = [yybackp cbackp ivbackp kp dp ap];

y = [gy gcc giv tby h yy c invs k1];
yp = [gyp gccp givp tbyp hp yyp cp invsp k1p];

%variables to substitute from levels to logs
logvar = [yyback cback ivback k d a gy gcc giv h yy c invs k1];
logvarp = [yybackp cbackp ivbackp kp dp ap gyp gccp givp hp yyp cp invsp k1p];

lv = [logvar logvarp];

% For the log-approximation: Make f a function of the logarithm of the state 
% and control vector
f = subs(f, lv, exp(lv));
logTransX = ones(1,6);
logTransY = ones(1,8);
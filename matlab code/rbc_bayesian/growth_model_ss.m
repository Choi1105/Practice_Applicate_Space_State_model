function [GAMA, DELTA, ALFA, PSSI, OMEGA, G, SIGMAG, RHOG, SIGMAA, RHOA, PHI, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, la, lap, a, ap, tby, tbyp, g, gp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, gback, gbackp, nx, ny, ne, sig, eta]=growth_model_ss(b);
%[GAMA, DELTA, ALFA, PSSI, OMEGA, G, SIGMAG, RHOG, SIGMAA, RHOA, PHI, RSTAR, DBAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, iv, ivp, tb, tbp, la, lap, a, ap, tby, tbyp, g, gp, gy, gyp, gc, gcp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, gback, gbackp]=growth_ss(b) produces the deep structural parameters and computes the steady state of a 
%small open economy described in the paper ``Real Business Cycles in Emerging Countries?,'' 
%by Javier Garcia-Cicco, Roberto Pancrazi, and Martin Uribe (forthcoming AER). 
%
%Input: vector b of estimated structural parameters
%
%(c)Martin Uribe
%
%September 2009

nx            = 8;      %Number of state variables
ny            = 9;      %Number of control variables
ne            = 2;      %Number of shocks
sig           = 1;

% PSSI = 0.001;%parameter governing the debt elasticity of the interest rate. 
% DBAR = 0.007; %Steady state of detrended external debt
GAMA = 2; %intertemporal elasticity of substitution

% G = 1.0107; %Gross growth rate of output
SIGMAG = b(1); %STD of innovation in permanent technology shock
RHOG = b(2); %Serial correlation of innovation in permanent technology shock
SIGMAA = b(3); %STD of innovation in transitory technology shock
RHOA = b(4); %Serial correlation of transitory technology shock
PHI = b(5); %Adjustment cost parameter/ 공장 증설 비용
G = b(6); %Adjustment cost parameter/ Growth rate
PSSI = b(7); %Adjustment cost parameter/ 부채에 대해 반응하는 이자율 정도

DELTA = 1.03^4-1;%Depreciation rate

ALFA = 0.32; %Capital elasticity of the production function

OMEGA = 1.6; %exponent of labor in utility function

THETA = 1.4 * OMEGA; %Preference parameter

RSTAR = 1.1; %parameter of the interest rate function
tby   = 0.003; %trade balance to output ratio 
BETTA = 1/RSTAR * G^GAMA; %World interest rate

% BETTA = 0.98^4;%0.98;%discount factor
% RSTAR = 1/BETTA * G^GAMA; %World interest rate

r=RSTAR;  %Country interest rate

% d =  DBAR; %foreign debt

k_over_gh = 	((G^GAMA/BETTA - 1 + DELTA) / ALFA)^(1/(ALFA-1)); %K/(G*H)

h = ((1-ALFA)* G * k_over_gh^ALFA / THETA)^(1/(OMEGA-1)); %hours
% h = ((1-ALFA)* G * k_over_gh^ALFA / 1)^(1/(OMEGA-1)); %hours

k = k_over_gh * G * h; %capital

invs = (G-1+DELTA) * k; %investment

yy = k^ALFA * (h*G)^(1-ALFA); %output

% YY = yy; %parameter of the premium function

% c = (G/r-1) * d + yy - invs; %Consumption

% tb = yy - c - invs; %Trade balance

% tby = tb / yy;
tb = tby * yy; %detrended trade balance

d = -tb/(G/r-1); %net external debt detrended

DY = d/yy; %parameter of interest-rate function, equal to the steady state level of net foreign debt. %ssg

c = (G/r-1) * d + yy - invs; %Consumption

la = (c - THETA/OMEGA*h^OMEGA)^(-GAMA); %marginal utility of wealth
% la = (c - 1/OMEGA*h^OMEGA)^(-GAMA); %marginal utility of wealth

k1 = k; %Auxiliary variable

a = 1; %productivity shock 

g = G; %Growth rate of nonstationary productivity shock

eta=[0 0 0 0 0 0 SIGMAG 0;  % Matrix defining driving force / shock : setting 가능
    0 0 0 0 0 0 0 SIGMAA]'; % sigmag = non stationary shock / sigmaa = stationary shock
                            

%Log variables
c = log(c);
k = log(k);
k1 = log(k1);
invs = log(invs);
h = log(h);
d = log(d);
la = log(la);
a = log(a);
g = log(g);
gy = g;
gcc = g;
giv = g;
yy = log(yy);
yyback = yy;
cback = c;
ivback = invs;
gback = g;

%Next-period variables
cp=c;
kp=k;
k1p=k;
invsp= invs;
hp=h;
dp=d; 
lap=la;
ap=a;
gp = g;
tbp = tb;
tbyp = tby;
gyp=gy;
gccp=gcc;
givp=giv;
yyp = yy;
yybackp=yy;
cbackp=cp;
ivbackp=invsp;
gbackp = g;
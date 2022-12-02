function [GAMA, DELTA, ALFA, PSSI, OMEGA, SIGMAA, RHOA, RSTAR, BETTA, THETA,...
    c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, YY, la, lap,...
    a, ap, tby, tbyp, gy, gyp, gcc, gccp, giv, givp, yy, yyp,...
    yyback, yybackp, cback, cbackp, ivback, ivbackp, nx, ny, ne, sig, eta]=growth_model_ss(b);
%Input: vector b of estimated structural parameters
%
%(c)Sanha Noh
%
%September 2021

nx            = 6;      %Number of state variables
ny            = 9;      %Number of control variables
ne            = 1;      %Number of shocks
sig           = 1;

% PSSI = 0.001;%parameter governing the debt elasticity of the interest rate. 
% DBAR = 0.007; %Steady state of detrended external debt
GAMA = 2; %intertemporal elasticity of substitution

% G = 1.0107; %Gross growth rate of output
% SIGMAG = 0; %STD of innovation in permanent technology shock
% RHOG   = 0; %Serial correlation of innovation in permanent technology shock
SIGMAA = b(1); %STD of innovation in transitory technology shock
RHOA   = b(2); %Serial correlation of transitory technology shock
% PHI    = 0; %Adjustment cost parameter
% G      = 1; %Adjustment cost parameter
PSSI   = 0.001; %Adjustment cost parameter

DELTA = 1.03^4-1;%Depreciation rate

ALFA = 0.32; %Capital elasticity of the production function

OMEGA = 1.6; %exponent of labor in utility function

THETA = 1.4 * OMEGA; %Preference parameter

RSTAR = 1.1; %parameter of the interest rate function
tby   = 0.003; %trade balance to output ratio 
BETTA = 1/RSTAR; %World interest rate

% BETTA = 0.98^4;%0.98;%discount factor
% RSTAR = 1/BETTA * G^GAMA; %World interest rate

r=RSTAR;  %Country interest rate

% d =  DBAR; %foreign debt

k_over_gh = 	((1/BETTA - 1 + DELTA) / ALFA)^(1/(ALFA-1)); %K/(G*H)

h = ((1-ALFA) * k_over_gh^ALFA / THETA)^(1/(OMEGA-1)); %hours
% h = ((1-ALFA)* G * k_over_gh^ALFA / 1)^(1/(OMEGA-1)); %hours

k = k_over_gh * h; %capital

invs = DELTA * k; %investment

yy = k^ALFA * h^(1-ALFA); %output

YY = yy; %parameter of the premium function

% c = (G/r-1) * d + yy - invs; %Consumption

% tb = yy - c - invs; %Trade balance

% tby = tb / yy;
tb = tby * yy; %detrended trade balance

d = -tb/(1/r-1); %net external debt detrended

DY = d/yy; %parameter of interest-rate function, equal to the steady state level of net foreign debt. %ssg

c = (1/r-1) * d + yy - invs; %Consumption

la = (c - THETA/OMEGA*h^OMEGA)^(-GAMA); %marginal utility of wealth

k1 = k; %Auxiliary variable

a = 1; %productivity shock 

eta=[0 0 0 0 0 SIGMAA]'; %Matrix defining driving force

%Log variables
c    = log(c);
k    = log(k);
k1   = log(k1);
invs = log(invs);
h    = log(h);
d    = log(d);
la   = log(la);
a    = log(a);
gy   = log(1);
gcc  = log(1);
giv  = log(1);
yy   = log(yy);
yyback = yy;
cback  = c;
ivback = invs;

%Next-period variables
cp=c;
kp=k;
k1p=k;
invsp= invs;
hp=h;
dp=d; 
lap=la;
ap=a;
tbp = tb;
tbyp = tby;
gyp=gy;
gccp=gcc;
givp=giv;
yyp = yy;
yybackp=yy;
cbackp=cp;
ivbackp=invsp;
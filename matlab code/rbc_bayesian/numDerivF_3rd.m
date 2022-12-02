function [fxxx,fxxxp,fxxy,fxxyp,fxxpxp,fxxpy,fxxpyp,fxyy,fxyyp,fxypyp,fxpxpxp,fxpxpy,fxpxpyp,fxpyy,fxpyyp,fxpypyp,fyyy,fyyyp,fyypyp,fypypyp,errorMes] = numDerivF_3rd(params)
%Initializing the flag for errors
errorMes = 0;
 
%The steady state of the model
[GAMA, DELTA, ALFA, PSSI, OMEGA, G, SIGMAG, RHOG, SIGMAA, RHOA, PHI, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, la, lap, a, ap, tby, tbyp, g, gp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, gback, gbackp, nx, ny, ne, sig, eta]=growth_model_ss(params);;
 
%Setting dimension for the matrices
n      = ny+nx;
%For fx:
fxxx   = zeros(n,nx,nx,nx);
fxxxp  = zeros(n,nx,nx,nx);
fxxy   = zeros(n,nx,nx,ny);
fxxyp  = zeros(n,nx,nx,ny);
fxxpxp = zeros(n,nx,nx,nx);
fxxpy  = zeros(n,nx,nx,ny);
fxxpyp = zeros(n,nx,nx,ny);
fxyy   = zeros(n,nx,ny,ny);
fxyyp  = zeros(n,nx,ny,ny);
fxypyp = zeros(n,nx,ny,ny);
%For fxp:
fxpxpxp= zeros(n,nx,nx,nx);
fxpxpy = zeros(n,nx,nx,ny);
fxpxpyp= zeros(n,nx,nx,ny);
fxpyy  = zeros(n,nx,ny,ny);
fxpyyp = zeros(n,nx,ny,ny);
fxpypyp= zeros(n,nx,ny,ny);
%For fy:
fyyy   = zeros(n,ny,ny,ny);
fyyyp  = zeros(n,ny,ny,ny);
fyypyp = zeros(n,ny,ny,ny);
%For fyp:
fypypyp = zeros(n,ny,ny,ny);
 
% START DISPLAYING fxxx
fxxx(11,1,1,1)=-exp(-yyback)*exp(gback)*exp(yy);
fxxx(11,4,1,1)=exp(-yyback)*exp(gback)*exp(yy);
fxxx(11,1,4,1)=exp(-yyback)*exp(gback)*exp(yy);
fxxx(11,4,4,1)=-exp(-yyback)*exp(gback)*exp(yy);
fxxx(12,2,2,2)=-exp(-cback)*exp(c)*exp(gback);
fxxx(12,4,2,2)=exp(-cback)*exp(c)*exp(gback);
fxxx(12,2,4,2)=exp(-cback)*exp(c)*exp(gback);
fxxx(12,4,4,2)=-exp(-cback)*exp(c)*exp(gback);
fxxx(13,3,3,3)=-exp(-ivback)*exp(gback)*exp(invs);
fxxx(13,4,3,3)=exp(-ivback)*exp(gback)*exp(invs);
fxxx(13,3,4,3)=exp(-ivback)*exp(gback)*exp(invs);
fxxx(13,4,4,3)=-exp(-ivback)*exp(gback)*exp(invs);
fxxx(11,1,1,4)=exp(-yyback)*exp(gback)*exp(yy);
fxxx(11,4,1,4)=-exp(-yyback)*exp(gback)*exp(yy);
fxxx(12,2,2,4)=exp(-cback)*exp(c)*exp(gback);
fxxx(12,4,2,4)=-exp(-cback)*exp(c)*exp(gback);
fxxx(13,3,3,4)=exp(-ivback)*exp(gback)*exp(invs);
fxxx(13,4,3,4)=-exp(-ivback)*exp(gback)*exp(invs);
fxxx(11,1,4,4)=-exp(-yyback)*exp(gback)*exp(yy);
fxxx(12,2,4,4)=-exp(-cback)*exp(c)*exp(gback);
fxxx(13,3,4,4)=-exp(-ivback)*exp(gback)*exp(invs);
fxxx(11,4,4,4)=exp(-yyback)*exp(gback)*exp(yy);
fxxx(12,4,4,4)=exp(-cback)*exp(c)*exp(gback);
fxxx(13,4,4,4)=exp(-ivback)*exp(gback)*exp(invs);
fxxx(1,5,5,5)=(PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxxx(2,5,5,5)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA) + 3*ALFA*exp(2*k)*exp(a)*exp(k)^(ALFA - 2)*(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA) + ALFA*exp(3*k)*exp(a)*exp(k)^(ALFA - 3)*(ALFA - 1)*(ALFA - 2)*(exp(g)*exp(h))^(1 - ALFA);
fxxx(3,5,5,5)=exp(k)*(DELTA - 1);
fxxx(5,5,5,5)=- 3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) - ALFA*exp(-3*h)*exp(3*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fxxx(6,5,5,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxx(10,5,5,5)=-exp(-yy)*((PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)));
fxxx(1,7,5,5)=PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxx(2,7,5,5)=- (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(2*k)*exp(a)*exp(g)*exp(h)*exp(k)^(ALFA - 2)*(ALFA - 1)^2)/(exp(g)*exp(h))^ALFA;
fxxx(5,7,5,5)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA + (ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^3)/exp(g)^ALFA;
fxxx(6,7,5,5)=-(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxx(10,7,5,5)=-PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxx(2,8,5,5)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA) + ALFA*exp(2*k)*exp(a)*exp(k)^(ALFA - 2)*(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxxx(5,8,5,5)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(1,5,7,5)=PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxx(2,5,7,5)=- (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(2*k)*exp(a)*exp(g)*exp(h)*exp(k)^(ALFA - 2)*(ALFA - 1)^2)/(exp(g)*exp(h))^ALFA;
fxxx(5,5,7,5)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA + (ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^3)/exp(g)^ALFA;
fxxx(6,5,7,5)=-(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxx(10,5,7,5)=-PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxx(1,7,7,5)=-2*PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxx(2,7,7,5)=(ALFA^2*exp(2*g)*exp(2*h)*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,7,7,5)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA^2*exp(2*g)*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^(ALFA + 1);
fxxx(6,7,7,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxx(10,7,7,5)=2*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxx(2,8,7,5)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,8,7,5)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,5,8,5)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA) + ALFA*exp(2*k)*exp(a)*exp(k)^(ALFA - 2)*(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxxx(5,5,8,5)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,7,8,5)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,7,8,5)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,8,8,5)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxxx(5,8,8,5)=-ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(1,6,6,6)=exp(d);
fxxx(1,5,5,7)=PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxx(2,5,5,7)=- (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(2*k)*exp(a)*exp(g)*exp(h)*exp(k)^(ALFA - 2)*(ALFA - 1)^2)/(exp(g)*exp(h))^ALFA;
fxxx(5,5,5,7)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA + (ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^3)/exp(g)^ALFA;
fxxx(6,5,5,7)=-(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxx(10,5,5,7)=-PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxx(1,7,5,7)=-2*PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxx(2,7,5,7)=(ALFA^2*exp(2*g)*exp(2*h)*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,7,5,7)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA^2*exp(2*g)*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^(ALFA + 1);
fxxx(6,7,5,7)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxx(10,7,5,7)=2*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxx(2,8,5,7)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,8,5,7)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(1,5,7,7)=-2*PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxx(2,5,7,7)=(ALFA^2*exp(2*g)*exp(2*h)*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,5,7,7)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA^2*exp(2*g)*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^(ALFA + 1);
fxxx(6,5,7,7)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxx(10,5,7,7)=2*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxx(1,7,7,7)=3*PHI*exp(2*g)*exp(-k)*exp(2*kp) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1));
fxxx(2,7,7,7)=(3*ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(3*g)*exp(3*h)*exp(a)*exp(k)^ALFA*(ALFA - 1)*(ALFA + 1))/(exp(g)*exp(h))^(ALFA + 2);
fxxx(3,7,7,7)=exp(g)*exp(kp);
fxxx(4,7,7,7)=(3*BETTA*GAMA*exp(2*g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*exp(g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*exp(3*g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^(GAMA + 3)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxx(5,7,7,7)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA - (3*ALFA*exp(2*g)*exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^(ALFA + 1) + (ALFA*exp(3*g)*exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2*(ALFA + 1))/exp(g)^(ALFA + 2);
fxxx(6,7,7,7)=(BETTA*GAMA*exp(g)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA - (3*BETTA*GAMA*exp(2*g)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*GAMA*exp(3*g)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 3)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxx(10,7,7,7)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - 3*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxx(17,7,7,7)=exp(g);
fxxx(2,8,7,7)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,8,7,7)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA*exp(2*g)*exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^(ALFA + 1);
fxxx(2,5,8,7)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,5,8,7)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,7,8,7)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,7,8,7)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA*exp(2*g)*exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^(ALFA + 1);
fxxx(2,8,8,7)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,8,8,7)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,5,5,8)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA) + ALFA*exp(2*k)*exp(a)*exp(k)^(ALFA - 2)*(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxxx(5,5,5,8)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,7,5,8)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,7,5,8)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,8,5,8)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxxx(5,8,5,8)=-ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,5,7,8)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,5,7,8)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,7,7,8)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,7,7,8)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA*exp(2*g)*exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^(ALFA + 1);
fxxx(2,8,7,8)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,8,7,8)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,5,8,8)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxxx(5,5,8,8)=-ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,7,8,8)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxx(5,7,8,8)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA;
fxxx(2,8,8,8)=exp(a)*exp(k)^ALFA*(exp(g)*exp(h))^(1 - ALFA);
fxxx(5,8,8,8)=-exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
% END DISPLAYING fxxx
if any(any(any(any(isnan(fxxx),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxx),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxx),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxxp
fxxxp(1,5,5,5)=PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxxp(6,5,5,5)=-(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxxp(10,5,5,5)=-PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxxp(1,7,5,5)=-2*PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxxp(6,7,5,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxxp(10,7,5,5)=2*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxxp(1,5,7,5)=-2*PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxxp(6,5,7,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxxp(10,5,7,5)=2*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxxp(1,7,7,5)=3*PHI*exp(2*g)*exp(-k)*exp(2*kp) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxxxp(3,7,7,5)=exp(g)*exp(kp);
fxxxp(6,7,7,5)=(ALFA*BETTA*GAMA*exp(2*g)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1)*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxxp(10,7,7,5)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - 3*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxxp(1,7,7,6)=(PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1));
fxxxp(4,7,7,6)=(BETTA*GAMA*PSSI*exp(2*g)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*PSSI*exp(-yy)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxxp(6,7,7,7)=(BETTA*GAMA*exp(2*g)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA)*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*exp(g)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxxp(6,7,7,8)=(ALFA*BETTA*GAMA*exp(2*g)*exp(ap)*(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (ALFA*BETTA*GAMA*exp(ap)*exp(g)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxxxp
if any(any(any(any(isnan(fxxxp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxxp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxxp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxy
fxxy(2,5,5,5)=- (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(2*k)*exp(a)*exp(g)*exp(h)*exp(k)^(ALFA - 2)*(ALFA - 1)^2)/(exp(g)*exp(h))^ALFA;
fxxy(5,5,5,5)=3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) + ALFA*exp(-3*h)*exp(3*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fxxy(6,5,5,5)=-(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(2,7,5,5)=(ALFA^2*exp(2*g)*exp(2*h)*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxy(5,7,5,5)=- (ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^3)/exp(g)^ALFA;
fxxy(6,7,5,5)=(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(2,8,5,5)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxy(5,8,5,5)=ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxy(2,5,7,5)=(ALFA^2*exp(2*g)*exp(2*h)*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxy(5,5,7,5)=- (ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^3)/exp(g)^ALFA;
fxxy(6,5,7,5)=(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(2,7,7,5)=(3*ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(3*g)*exp(3*h)*exp(a)*exp(k)^ALFA*(ALFA - 1)*(ALFA + 1))/(exp(g)*exp(h))^(ALFA + 2);
fxxy(5,7,7,5)=(ALFA^2*exp(2*g)*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^(ALFA + 1) - (ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxy(6,7,7,5)=-(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(2,8,7,5)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxy(5,8,7,5)=-(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxy(2,5,8,5)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxy(5,5,8,5)=ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxy(2,7,8,5)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxy(5,7,8,5)=-(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxxy(2,8,8,5)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxxy(5,8,8,5)=ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxy(11,1,1,6)=exp(-yyback)*exp(gback)*exp(yy);
fxxy(11,4,1,6)=-exp(-yyback)*exp(gback)*exp(yy);
fxxy(11,1,4,6)=-exp(-yyback)*exp(gback)*exp(yy);
fxxy(11,4,4,6)=exp(-yyback)*exp(gback)*exp(yy);
fxxy(10,5,5,6)=exp(-yy)*((PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) + PHI*exp(2*g)*exp(-k)*exp(2*kp));
fxxy(10,7,5,6)=-PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxy(10,5,7,6)=-PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxy(1,7,7,6)=-(PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fxxy(4,7,7,6)=(BETTA*GAMA*PSSI*exp(-yy)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*PSSI*exp(2*g)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxy(10,7,7,6)=PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy) - PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxxy(12,2,2,7)=exp(-cback)*exp(c)*exp(gback);
fxxy(12,4,2,7)=-exp(-cback)*exp(c)*exp(gback);
fxxy(12,2,4,7)=-exp(-cback)*exp(c)*exp(gback);
fxxy(12,4,4,7)=exp(-cback)*exp(c)*exp(gback);
fxxy(6,5,5,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(6,7,5,7)=-(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(6,5,7,7)=-(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(6,7,7,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxy(13,3,3,8)=exp(-ivback)*exp(gback)*exp(invs);
fxxy(13,4,3,8)=-exp(-ivback)*exp(gback)*exp(invs);
fxxy(13,3,4,8)=-exp(-ivback)*exp(gback)*exp(invs);
fxxy(13,4,4,8)=exp(-ivback)*exp(gback)*exp(invs);
fxxy(6,7,7,9)=(BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*PHI*exp(2*g)*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxxy
if any(any(any(any(isnan(fxxy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxyp
fxxyp(4,7,7,5)=(BETTA*GAMA^2*THETA*exp(2*g)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA^2*THETA*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxyp(6,7,7,5)=(BETTA*GAMA^2*THETA*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA^2*THETA*exp(2*g)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*exp(2*g)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1)*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxyp(4,7,7,7)=(BETTA*GAMA^2*exp(cp)*exp(g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA^2*exp(2*g)*exp(cp)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxyp(6,7,7,7)=(BETTA*GAMA^2*exp(2*g)*exp(cp)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA^2*exp(cp)*exp(g)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxyp(6,7,7,9)=(BETTA*GAMA*PHI*exp(2*g)*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxxyp
if any(any(any(any(isnan(fxxyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxpxp
fxxpxp(1,5,5,5)=-2*PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxpxp(6,5,5,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxpxp(10,5,5,5)=2*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxpxp(1,7,5,5)=3*PHI*exp(2*g)*exp(-k)*exp(2*kp) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxxpxp(3,7,5,5)=exp(g)*exp(kp);
fxxpxp(6,7,5,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(g)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxpxp(10,7,5,5)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - 3*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxpxp(6,7,7,5)=(BETTA*GAMA*exp(g)*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxpxp(6,7,8,5)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxpxp(1,7,6,6)=(3*PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (2*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)) + (PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fxxpxp(4,7,6,6)=- (BETTA*GAMA*PSSI*exp(-yy)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*PSSI*exp(2*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxpxp(6,7,5,7)=(ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(g)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxpxp(6,7,7,7)=-(BETTA*GAMA*exp(g)*(2*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) + (ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxpxp(6,7,8,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxpxp(6,7,5,8)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxpxp(6,7,7,8)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxpxp(6,7,8,8)=-(ALFA*BETTA*GAMA*exp(ap)*exp(g)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxxpxp
if any(any(any(any(isnan(fxxpxp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxpxp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxpxp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxpy
fxxpy(6,5,5,5)=(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxpy(6,7,5,5)=-(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxpy(10,5,5,6)=-PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxpy(10,7,5,6)=PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy) - PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxxpy(1,7,6,6)=(2*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (2*PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fxxpy(4,7,6,6)=(BETTA*GAMA*PSSI*exp(-yy)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*GAMA*PSSI*exp(2*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxpy(6,5,5,7)=-(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxpy(6,7,5,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxxpy(6,7,7,9)=(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxxpy
if any(any(any(any(isnan(fxxpy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxpy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxpy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxpyp
fxxpyp(6,7,5,5)=(ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(g)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA^2*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxpyp(4,7,6,5)=-(BETTA*GAMA^2*PSSI*THETA*exp(-yy)*exp(dp)*exp(g)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxpyp(6,7,7,5)=- (BETTA*GAMA*exp(g)*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA^2*THETA*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxpyp(6,7,8,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA^2*THETA*exp(ap)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxpyp(6,7,5,7)=(ALFA*BETTA*GAMA^2*exp(-kp)*exp(ap)*exp(cp)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxpyp(4,7,6,7)=(BETTA*GAMA^2*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxpyp(6,7,7,7)=(BETTA*GAMA^2*exp(cp)*exp(g)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxpyp(6,7,8,7)=(ALFA*BETTA*GAMA^2*exp(ap)*exp(cp)*exp(g)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxxpyp(6,7,7,9)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxxpyp
if any(any(any(any(isnan(fxxpyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxpyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxpyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxyy
fxyy(2,5,5,5)=(ALFA^2*exp(2*g)*exp(2*h)*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxyy(5,5,5,5)=- 3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) - ALFA*exp(-3*h)*exp(3*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fxyy(6,5,5,5)=(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) + (GAMA*PHI*THETA*exp(2*h)*exp(-k)*exp(g)*exp(kp)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) + (GAMA*PHI*THETA^2*exp(2*h)*exp(-k)*exp(g)*exp(kp)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxyy(2,7,5,5)=(3*ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(3*g)*exp(3*h)*exp(a)*exp(k)^ALFA*(ALFA - 1)*(ALFA + 1))/(exp(g)*exp(h))^(ALFA + 2);
fxyy(5,7,5,5)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA + (ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^3)/exp(g)^ALFA;
fxyy(6,7,5,5)=- (GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*PHI*THETA*exp(2*h)*exp(-k)*exp(g)*exp(kp)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*PHI*THETA^2*exp(2*h)*exp(-k)*exp(g)*exp(kp)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxyy(2,8,5,5)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxyy(5,8,5,5)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxyy(6,5,7,5)=-(GAMA*PHI*THETA*exp(-k)*exp(c)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxyy(6,7,7,5)=(GAMA*PHI*THETA*exp(-k)*exp(c)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxyy(11,1,6,6)=-exp(-yyback)*exp(gback)*exp(yy);
fxyy(11,4,6,6)=exp(-yyback)*exp(gback)*exp(yy);
fxyy(10,5,6,6)=-exp(-yy)*((PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)));
fxyy(1,7,6,6)=(PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (2*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 + (PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fxyy(4,7,6,6)=- (BETTA*GAMA*PSSI*exp(-yy)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*PSSI*exp(2*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxyy(10,7,6,6)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxyy(6,5,5,7)=-(GAMA*PHI*THETA*exp(-k)*exp(c)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxyy(6,7,5,7)=(GAMA*PHI*THETA*exp(-k)*exp(c)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxyy(12,2,7,7)=-exp(-cback)*exp(c)*exp(gback);
fxyy(12,4,7,7)=exp(-cback)*exp(c)*exp(gback);
fxyy(6,5,7,7)=(GAMA*PHI*exp(2*c)*exp(-k)*exp(g)*exp(kp)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxyy(6,7,7,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*PHI*exp(2*c)*exp(-k)*exp(g)*exp(kp)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxyy(13,3,8,8)=-exp(-ivback)*exp(gback)*exp(invs);
fxyy(13,4,8,8)=exp(-ivback)*exp(gback)*exp(invs);
fxyy(6,7,9,9)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxyy
if any(any(any(any(isnan(fxyy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxyy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxyy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxyyp
fxyyp(4,7,6,5)=(BETTA*GAMA^2*PSSI*THETA*exp(-yy)*exp(dp)*exp(g)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxyyp(6,7,9,5)=(BETTA*GAMA^2*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxyyp(4,7,6,7)=-(BETTA*GAMA^2*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxyyp(6,7,9,7)=-(BETTA*GAMA^2*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxyyp(6,7,9,9)=(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxyyp
if any(any(any(any(isnan(fxyyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxyyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxyyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxypyp
fxypyp(4,7,5,5)=- (BETTA*GAMA^2*THETA*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA^2*THETA^2*exp(2*hp)*exp(g)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA^2*THETA*exp(2*hp)*exp(g)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(OMEGA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxypyp(6,7,5,5)=(BETTA*GAMA^2*THETA*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA^2*THETA^2*exp(2*hp)*exp(g)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA^2*THETA*exp(2*hp)*exp(g)*exp(hp)^(OMEGA - 2)*(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(g)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (2*ALFA*BETTA*GAMA^2*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxypyp(4,7,7,5)=(BETTA*GAMA^2*THETA*exp(cp)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxypyp(6,7,7,5)=- (BETTA*GAMA^2*THETA*exp(cp)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (ALFA*BETTA*GAMA^2*exp(-kp)*exp(ap)*exp(cp)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxypyp(6,7,9,5)=-(BETTA*GAMA^2*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxypyp(4,7,5,7)=(BETTA*GAMA^2*THETA*exp(cp)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxypyp(6,7,5,7)=- (BETTA*GAMA^2*THETA*exp(cp)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (ALFA*BETTA*GAMA^2*exp(-kp)*exp(ap)*exp(cp)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxypyp(4,7,7,7)=(BETTA*GAMA^2*exp(cp)*exp(g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA^2*exp(2*cp)*exp(g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxypyp(6,7,7,7)=(BETTA*GAMA^2*exp(2*cp)*exp(g)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA^2*exp(cp)*exp(g)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxypyp(6,7,9,7)=(BETTA*GAMA^2*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxypyp(6,7,5,9)=-(BETTA*GAMA^2*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxypyp(6,7,7,9)=(BETTA*GAMA^2*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxypyp(6,7,9,9)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxypyp
if any(any(any(any(isnan(fxypyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxypyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxypyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxpxp
fxpxpxp(14,1,1,1)=-exp(yybackp);
fxpxpxp(15,2,2,2)=-exp(cbackp);
fxpxpxp(16,3,3,3)=-exp(ivbackp);
fxpxpxp(17,4,4,4)=-exp(gbackp);
fxpxpxp(1,5,5,5)=3*PHI*exp(2*g)*exp(-k)*exp(2*kp) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxpxpxp(3,5,5,5)=exp(g)*exp(kp);
fxpxpxp(6,5,5,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA - (3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2));
fxpxpxp(7,5,5,5)=exp(kp);
fxpxpxp(10,5,5,5)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - 3*PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxpxpxp(6,7,5,5)=-(BETTA*((ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA - (3*ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) + (ALFA^2*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,8,5,5)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,5,7,5)=(3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2));
fxpxpxp(6,7,7,5)=(BETTA*((ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA - (3*ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) + (ALFA^2*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,8,7,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpxp(6,5,8,5)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,7,8,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpxp(6,8,8,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(1,6,6,6)=(7*PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (12*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (6*PSSI^2*exp(4*dp)*exp(-3*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)) + (6*PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 + (PSSI*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 + (6*PSSI^3*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^4;
fxpxpxp(4,6,6,6)=(3*BETTA*PSSI*exp(2*dp)*exp(-2*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(3*dp)*exp(-3*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,5,5,7)=(3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2));
fxpxpxp(6,7,5,7)=(BETTA*((ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA - (3*ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) + (ALFA^2*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,8,5,7)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpxp(6,5,7,7)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2));
fxpxpxp(6,7,7,7)=(BETTA*(4*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) + (3*ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA - (ALFA^2*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,8,7,7)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,5,8,7)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpxp(6,7,8,7)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,8,8,7)=-(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,5,5,8)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,7,5,8)=-(BETTA*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,8,5,8)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,5,7,8)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpxp(6,7,7,8)=(BETTA*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,8,7,8)=-(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,5,8,8)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,7,8,8)=-(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpxp(6,8,8,8)=(ALFA*BETTA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpxpxp
if any(any(any(any(isnan(fxpxpxp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpxpxp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpxpxp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxpy
fxpxpy(6,5,5,5)=-(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxpxpy(10,5,5,6)=PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy) - PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxpxpy(1,6,6,6)=(10*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 + (6*PSSI^2*exp(4*dp)*exp(-3*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (4*PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (5*PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (PSSI*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (6*PSSI^3*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^4;
fxpxpy(4,6,6,6)=- (3*BETTA*PSSI*exp(2*dp)*exp(-2*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*PSSI*exp(3*dp)*exp(-3*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpy(6,5,5,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxpxpy(6,7,7,9)=-(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpxpy
if any(any(any(any(isnan(fxpxpy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpxpy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpxpy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxpyp
fxpxpyp(6,5,5,5)=(3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (ALFA^2*BETTA*GAMA*THETA*exp(2*gp)*exp(3*hp)*exp(-2*kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpyp(6,7,5,5)=(BETTA*((ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA - (3*ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) + (ALFA^2*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,8,5,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(4,6,6,5)=(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PSSI*THETA*exp(2*dp)*exp(-2*yy)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,5,7,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*GAMA*THETA*exp(2*gp)*exp(3*hp)*exp(-2*kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpyp(6,7,7,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(2*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) + (ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*((ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA - (3*ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) + (ALFA^2*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpyp(6,8,7,5)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,5,8,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,7,8,5)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,8,8,5)=(ALFA*BETTA*GAMA*THETA*exp(ap)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,5,5,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxpyp(6,7,5,7)=(BETTA*GAMA*exp(cp)*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,8,5,7)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(4,6,6,7)=- (BETTA*GAMA*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PSSI*exp(2*dp)*exp(-2*yy)*exp(cp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,5,7,7)=(ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,7,7,7)=-(BETTA*GAMA*exp(cp)*(2*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) + (ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,8,7,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,5,8,7)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,7,8,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxpyp(6,8,8,7)=-(ALFA*BETTA*GAMA*exp(ap)*exp(cp)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,7,7,9)=(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpxpyp
if any(any(any(any(isnan(fxpxpyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpxpyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpxpyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpyy
fxpyy(6,5,5,5)=- (GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*PHI*THETA*exp(2*h)*exp(-k)*exp(g)*exp(kp)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*PHI*THETA^2*exp(2*h)*exp(-k)*exp(g)*exp(kp)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxpyy(6,5,7,5)=(GAMA*PHI*THETA*exp(-k)*exp(c)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxpyy(10,5,6,6)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxpyy(1,6,6,6)=(2*PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (6*PSSI^2*exp(4*dp)*exp(-3*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (8*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 + (4*PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 + (PSSI*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 + (6*PSSI^3*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^4;
fxpyy(4,6,6,6)=(3*BETTA*PSSI*exp(2*dp)*exp(-2*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(3*dp)*exp(-3*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpyy(6,5,5,7)=(GAMA*PHI*THETA*exp(-k)*exp(c)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxpyy(6,5,7,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*PHI*exp(2*c)*exp(-k)*exp(g)*exp(kp)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fxpyy(6,7,9,9)=(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpyy
if any(any(any(any(isnan(fxpyy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpyy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpyy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpyyp
fxpyyp(4,6,6,5)=- (BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PSSI*THETA*exp(2*dp)*exp(-2*yy)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyyp(6,7,9,5)=-(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyyp(4,6,6,7)=(BETTA*GAMA*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PSSI*exp(2*dp)*exp(-2*yy)*exp(cp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyyp(6,7,9,7)=(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyyp(6,7,9,9)=-(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpyyp
if any(any(any(any(isnan(fxpyyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpyyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpyyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpypyp
fxpypyp(6,5,5,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)) + (3*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (2*ALFA^2*BETTA*GAMA*THETA*exp(2*gp)*exp(3*hp)*exp(-2*kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*THETA*exp(3*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 2)*(ALFA - 1)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (ALFA*BETTA*GAMA*THETA^2*exp(3*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(2*OMEGA - 2)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpypyp(4,6,5,5)=(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PSSI*THETA^2*exp(2*hp)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*PSSI*THETA*exp(2*hp)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 2)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,7,5,5)=(2*BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*((ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA - (3*ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) + (ALFA^2*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,8,5,5)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (ALFA*BETTA*GAMA*THETA*exp(ap)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (ALFA*BETTA*GAMA*THETA^2*exp(2*hp)*exp(ap)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(ap)*exp(hp)^(OMEGA - 2)*(OMEGA - 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpypyp(6,5,7,5)=(ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpypyp(4,6,7,5)=-(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(cp)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,7,7,5)=- (BETTA*GAMA*exp(cp)*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,8,7,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(ap)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,7,9,5)=(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,5,5,7)=(ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpypyp(4,6,5,7)=-(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(cp)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,7,5,7)=- (BETTA*GAMA*exp(cp)*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,8,5,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(ap)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,5,7,7)=(ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpypyp(4,6,7,7)=(BETTA*GAMA*PSSI*exp(2*cp)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,7,7,7)=(BETTA*GAMA*exp(2*cp)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*exp(cp)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,8,7,7)=(ALFA*BETTA*GAMA*exp(2*cp)*exp(ap)*(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (ALFA*BETTA*GAMA*exp(ap)*exp(cp)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,7,9,7)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,7,5,9)=(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,7,7,9)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,7,9,9)=(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpypyp
if any(any(any(any(isnan(fxpypyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpypyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpypyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyyy
fyyy(11,1,1,1)=-exp(gy);
fyyy(12,2,2,2)=-exp(gcc);
fyyy(13,3,3,3)=-exp(giv);
fyyy(2,5,5,5)=(3*ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA - (ALFA*exp(3*g)*exp(3*h)*exp(a)*exp(k)^ALFA*(ALFA - 1)*(ALFA + 1))/(exp(g)*exp(h))^(ALFA + 2);
fyyy(4,5,5,5)=- (GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(3*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA^2*exp(3*h)*exp(h)^(2*OMEGA - 3)*(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^2*exp(3*h)*exp(h)^(OMEGA - 1)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^3*exp(3*h)*exp(h)^(2*OMEGA - 2)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(5,5,5,5)=3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - 3*THETA*exp(2*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2) - THETA*exp(3*h)*exp(h)^(OMEGA - 4)*(OMEGA - 1)*(OMEGA - 2)*(OMEGA - 3) - THETA*exp(h)*exp(h)^(OMEGA - 2)*(OMEGA - 1) + ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) + ALFA*exp(-3*h)*exp(3*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fyyy(6,5,5,5)=(GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) + (3*GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) + (3*GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA*exp(3*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) + (GAMA*THETA^2*exp(3*h)*exp(h)^(2*OMEGA - 3)*(2*OMEGA - 2)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^3*exp(3*h)*exp(h)^(2*OMEGA - 2)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA^2*exp(3*h)*exp(h)^(OMEGA - 1)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(4,7,5,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(6,7,5,5)=- (GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) - (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(4,5,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(6,5,7,5)=- (GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) - (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(4,7,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,7,7,5)=(GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) - (GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(1,6,6,6)=(6*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - exp(yy) + (6*PSSI^2*exp(4*dp)*exp(-3*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (3*PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (PSSI*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (6*PSSI^3*exp(4*dp)*exp(-3*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^4;
fyyy(2,6,6,6)=-exp(yy);
fyyy(4,6,6,6)=- (3*BETTA*PSSI*exp(2*dp)*exp(-2*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*PSSI*exp(3*dp)*exp(-3*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fyyy(10,6,6,6)=exp(-yy)*(exp(c) + exp(invs) - exp(yy) + (PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2) + 1;
fyyy(11,6,6,6)=exp(-yyback)*exp(gback)*exp(yy);
fyyy(14,6,6,6)=exp(yy);
fyyy(10,7,6,6)=-exp(-yy)*exp(c);
fyyy(10,8,6,6)=-exp(-yy)*exp(invs);
fyyy(10,6,7,6)=-exp(-yy)*exp(c);
fyyy(10,7,7,6)=exp(-yy)*exp(c);
fyyy(10,6,8,6)=-exp(-yy)*exp(invs);
fyyy(10,8,8,6)=exp(-yy)*exp(invs);
fyyy(4,5,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(6,5,5,7)=- (GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) - (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(4,7,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,7,5,7)=(GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) - (GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(10,6,6,7)=-exp(-yy)*exp(c);
fyyy(10,7,6,7)=exp(-yy)*exp(c);
fyyy(4,5,7,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,5,7,7)=(GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) - (GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(10,6,7,7)=exp(-yy)*exp(c);
fyyy(1,7,7,7)=exp(c);
fyyy(4,7,7,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*exp(2*c)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*exp(3*c)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,7,7,7)=(3*GAMA*exp(2*c)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*exp(c)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*exp(3*c)*(GAMA + 1)*(GAMA + 2)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(10,7,7,7)=-exp(-yy)*exp(c);
fyyy(12,7,7,7)=exp(-cback)*exp(c)*exp(gback);
fyyy(15,7,7,7)=exp(c);
fyyy(10,6,6,8)=-exp(-yy)*exp(invs);
fyyy(10,8,6,8)=exp(-yy)*exp(invs);
fyyy(10,6,8,8)=exp(-yy)*exp(invs);
fyyy(1,8,8,8)=exp(invs);
fyyy(3,8,8,8)=-exp(invs);
fyyy(10,8,8,8)=-exp(-yy)*exp(invs);
fyyy(13,8,8,8)=exp(-ivback)*exp(gback)*exp(invs);
fyyy(16,8,8,8)=exp(invs);
fyyy(6,9,9,9)=-(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fyyy(7,9,9,9)=-exp(k1);
% END DISPLAYING fyyy
if any(any(any(any(isnan(fyyy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fyyy),1),2),3),4) ~= 0 || any(any(any(any(imag(fyyy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyyyp
fyyyp(4,6,6,5)=(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PSSI*THETA*exp(2*dp)*exp(-2*yy)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyyp(6,9,9,5)=(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyyp(4,6,6,7)=- (BETTA*GAMA*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PSSI*exp(2*dp)*exp(-2*yy)*exp(cp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyyp(6,9,9,7)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyyp(6,9,9,9)=(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fyyyp
if any(any(any(any(isnan(fyyyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fyyyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fyyyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyypyp
fyypyp(4,6,5,5)=- (BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PSSI*THETA^2*exp(2*hp)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*PSSI*THETA*exp(2*hp)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 2)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyypyp(6,9,5,5)=- (BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(OMEGA - 2)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PHI*THETA^2*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fyypyp(4,6,7,5)=(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(cp)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fyypyp(6,9,7,5)=(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fyypyp(6,9,9,5)=-(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyypyp(4,6,5,7)=(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(cp)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fyypyp(6,9,5,7)=(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fyypyp(4,6,7,7)=(BETTA*GAMA*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PSSI*exp(2*cp)*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fyypyp(6,9,7,7)=(BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PHI*exp(2*cp)*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fyypyp(6,9,9,7)=(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyypyp(6,9,5,9)=-(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyypyp(6,9,7,9)=(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyypyp(6,9,9,9)=-(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fyypyp
if any(any(any(any(isnan(fyypyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fyypyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fyypyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fypypyp
fypypyp(4,5,5,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (3*BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (3*BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*THETA*exp(3*hp)*exp(hp)^(OMEGA - 3)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(OMEGA - 1)*(OMEGA - 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(2*OMEGA - 3)*(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(OMEGA - 1)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*THETA^3*exp(3*hp)*exp(hp)^(2*OMEGA - 2)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3));
fypypyp(6,5,5,5)=(3*ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (3*BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (3*BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*THETA*exp(3*hp)*exp(hp)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(2*OMEGA - 3)*(2*OMEGA - 2)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA^3*exp(3*hp)*exp(hp)^(2*OMEGA - 2)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) - (ALFA^2*BETTA*exp(3*gp)*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 2)) - (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(OMEGA - 1)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (6*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (3*ALFA^2*BETTA*GAMA*THETA*exp(2*gp)*exp(3*hp)*exp(-2*kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (3*ALFA*BETTA*GAMA*THETA*exp(3*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 2)*(ALFA - 1)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (3*ALFA*BETTA*GAMA*THETA^2*exp(3*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(2*OMEGA - 2)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypypyp(4,7,5,5)=- (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3));
fypypyp(6,7,5,5)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) + (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypypyp(6,9,5,5)=(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(OMEGA - 2)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PHI*THETA^2*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(4,5,7,5)=- (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3));
fypypyp(6,5,7,5)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) + (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypypyp(4,7,7,5)=(BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,7,7,5)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypypyp(6,9,7,5)=-(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,5,9,5)=(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(OMEGA - 2)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PHI*THETA^2*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,7,9,5)=-(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,9,9,5)=(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(4,5,5,7)=- (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3));
fypypyp(6,5,5,7)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) + (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (ALFA^2*BETTA*GAMA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) + (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypypyp(4,7,5,7)=(BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,7,5,7)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypypyp(6,9,5,7)=-(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(4,5,7,7)=(BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,5,7,7)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypypyp(4,7,7,7)=(3*BETTA*GAMA*exp(2*cp)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*exp(cp)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*exp(3*cp)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3));
fypypyp(6,7,7,7)=(BETTA*GAMA*exp(cp)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (3*BETTA*GAMA*exp(2*cp)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*exp(3*cp)*(GAMA + 1)*(GAMA + 2)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3));
fypypyp(6,9,7,7)=(BETTA*GAMA*PHI*exp(2*cp)*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,5,9,7)=-(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,7,9,7)=(BETTA*GAMA*PHI*exp(2*cp)*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,9,9,7)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,5,5,9)=(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(OMEGA - 2)*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PHI*THETA^2*exp(2*gp)*exp(2*hp)*exp(-2*k1)*exp(2*k1p)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,7,5,9)=-(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,9,5,9)=(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,5,7,9)=-(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypypyp(6,7,7,9)=(BETTA*GAMA*PHI*exp(2*cp)*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,9,7,9)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,5,9,9)=(2*BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,7,9,9)=-(2*BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypypyp(6,9,9,9)=(4*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fypypyp
if any(any(any(any(isnan(fypypyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fypypyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fypypyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
if errorMes == 1
    fxxx   = NaN;
    fxxxp  = NaN;
    fxxy   = NaN;
    fxxyp  = NaN;
    fxxpxp = NaN;
    fxxpy  = NaN;
    fxxpyp = NaN;
    fxyy   = NaN;
    fxyyp  = NaN;
    fxypyp = NaN;
    fxpxpxp= NaN;
    fxpxpy = NaN;
    fxpxpyp= NaN;
    fxpyy  = NaN;
    fxpyyp = NaN;
    fxpypyp= NaN;
    fyyy   = NaN;
    fyyyp  = NaN;
    fyypyp = NaN;
    fypypyp= NaN;
end
end

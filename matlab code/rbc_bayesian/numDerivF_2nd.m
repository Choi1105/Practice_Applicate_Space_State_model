function [fxx,fxxp,fxy,fxyp,fxpx,fxpxp,fxpy,fxpyp,fyx,fyxp,fyy,fyyp,fypx,fypxp,fypy,fypyp,errorMes] = numDerivF_2nd(params)
 
%Initializing the flag for errors
errorMes = 0;
 
%The steady state of the model
[GAMA, DELTA, ALFA, PSSI, OMEGA, G, SIGMAG, RHOG, SIGMAA, RHOA, PHI, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, la, lap, a, ap, tby, tbyp, g, gp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, gback, gbackp, nx, ny, ne, sig, eta]=growth_model_ss(params);;
 
%Setting dimension for the matrices
n     = ny+nx;
fxx   = zeros(n,nx,nx);
fxxp  = zeros(n,nx,nx);
fxy   = zeros(n,nx,ny);
fxyp  = zeros(n,nx,ny);
fxpxp = zeros(n,nx,nx);
fxpy  = zeros(n,nx,ny);
fxpyp = zeros(n,nx,ny);
fyy   = zeros(n,ny,ny);
fyyp  = zeros(n,ny,ny);
fypyp = zeros(n,ny,ny);
 
% START DISPLAYING fxx
fxx(11,1,1)=exp(-yyback)*exp(gback)*exp(yy);
fxx(11,4,1)=-exp(-yyback)*exp(gback)*exp(yy);
fxx(12,2,2)=exp(-cback)*exp(c)*exp(gback);
fxx(12,4,2)=-exp(-cback)*exp(c)*exp(gback);
fxx(13,3,3)=exp(-ivback)*exp(gback)*exp(invs);
fxx(13,4,3)=-exp(-ivback)*exp(gback)*exp(invs);
fxx(11,1,4)=-exp(-yyback)*exp(gback)*exp(yy);
fxx(12,2,4)=-exp(-cback)*exp(c)*exp(gback);
fxx(13,3,4)=-exp(-ivback)*exp(gback)*exp(invs);
fxx(11,4,4)=exp(-yyback)*exp(gback)*exp(yy);
fxx(12,4,4)=exp(-cback)*exp(c)*exp(gback);
fxx(13,4,4)=exp(-ivback)*exp(gback)*exp(invs);
fxx(1,5,5)=(PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) + PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxx(2,5,5)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA) + ALFA*exp(2*k)*exp(a)*exp(k)^(ALFA - 2)*(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxx(3,5,5)=exp(k)*(DELTA - 1);
fxx(5,5,5)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxx(6,5,5)=-(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxx(10,5,5)=-exp(-yy)*((PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) + PHI*exp(2*g)*exp(-k)*exp(2*kp));
fxx(1,7,5)=-PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxx(2,7,5)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxx(5,7,5)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxx(6,7,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxx(10,7,5)=PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxx(2,8,5)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxx(5,8,5)=-ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxx(1,6,6)=exp(d);
fxx(1,5,7)=-PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxx(2,5,7)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxx(5,5,7)=(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxx(6,5,7)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxx(10,5,7)=PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxx(1,7,7)=PHI*exp(2*g)*exp(-k)*exp(2*kp) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1));
fxx(2,7,7)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxx(3,7,7)=exp(g)*exp(kp);
fxx(4,7,7)=(BETTA*GAMA*exp(2*g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*GAMA*exp(g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxx(5,7,7)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA - (ALFA*exp(2*g)*exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^(ALFA + 1);
fxx(6,7,7)=(BETTA*GAMA*exp(g)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA - (BETTA*GAMA*exp(2*g)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 2)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxx(10,7,7)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxx(17,7,7)=exp(g);
fxx(2,8,7)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxx(5,8,7)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA;
fxx(2,5,8)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fxx(5,5,8)=-ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxx(2,7,8)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxx(5,7,8)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA;
fxx(2,8,8)=exp(a)*exp(k)^ALFA*(exp(g)*exp(h))^(1 - ALFA);
fxx(5,8,8)=-exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
% END DISPLAYING fxx
if any(any(any(isnan(fxx),1),2),3) ~= 0 || any(any(any(isinf(fxx),1),2),3) ~= 0 || any(any(any(imag(fxx),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxp
fxxp(1,5,5)=-PHI*exp(2*g)*exp(-k)*exp(2*kp);
fxxp(6,5,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxxp(10,5,5)=PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxp(1,7,5)=PHI*exp(2*g)*exp(-k)*exp(2*kp) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxxp(3,7,5)=exp(g)*exp(kp);
fxxp(6,7,5)=- (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxxp(10,7,5)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxxp(1,7,6)=(PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1));
fxxp(4,7,6)=-(BETTA*GAMA*PSSI*exp(-yy)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxp(6,7,7)=-(BETTA*GAMA*exp(g)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxxp(6,7,8)=-(ALFA*BETTA*GAMA*exp(ap)*exp(g)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxxp
if any(any(any(isnan(fxxp),1),2),3) ~= 0 || any(any(any(isinf(fxxp),1),2),3) ~= 0 || any(any(any(imag(fxxp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxy
fxy(2,5,5)=-(ALFA*exp(a)*exp(g)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxy(5,5,5)=ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxy(6,5,5)=(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxy(2,7,5)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxy(5,7,5)=-(ALFA*exp(-h)*exp(a)*exp(g)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1)^2)/exp(g)^ALFA;
fxy(6,7,5)=-(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxy(2,8,5)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fxy(5,8,5)=ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxy(11,1,6)=-exp(-yyback)*exp(gback)*exp(yy);
fxy(11,4,6)=exp(-yyback)*exp(gback)*exp(yy);
fxy(10,5,6)=exp(-yy)*((PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)));
fxy(1,7,6)=-(PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fxy(4,7,6)=(BETTA*GAMA*PSSI*exp(-yy)*exp(dp)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxy(10,7,6)=-PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxy(12,2,7)=-exp(-cback)*exp(c)*exp(gback);
fxy(12,4,7)=exp(-cback)*exp(c)*exp(gback);
fxy(6,5,7)=-(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxy(6,7,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxy(13,3,8)=-exp(-ivback)*exp(gback)*exp(invs);
fxy(13,4,8)=exp(-ivback)*exp(gback)*exp(invs);
fxy(6,7,9)=(BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxy
if any(any(any(isnan(fxy),1),2),3) ~= 0 || any(any(any(isinf(fxy),1),2),3) ~= 0 || any(any(any(imag(fxy),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxyp
fxyp(4,7,5)=-(BETTA*GAMA^2*THETA*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxyp(6,7,5)=(BETTA*GAMA^2*THETA*exp(g)*exp(hp)*exp(hp)^(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(g)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxyp(4,7,7)=(BETTA*GAMA^2*exp(cp)*exp(g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxyp(6,7,7)=-(BETTA*GAMA^2*exp(cp)*exp(g)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxyp(6,7,9)=-(BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(g))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxyp
if any(any(any(isnan(fxyp),1),2),3) ~= 0 || any(any(any(isinf(fxyp),1),2),3) ~= 0 || any(any(any(imag(fxyp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxp
fxpxp(14,1,1)=-exp(yybackp);
fxpxp(15,2,2)=-exp(cbackp);
fxpxp(16,3,3)=-exp(ivbackp);
fxpxp(17,4,4)=-exp(gbackp);
fxpxp(1,5,5)=PHI*exp(2*g)*exp(-k)*exp(2*kp) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxpxp(3,5,5)=exp(g)*exp(kp);
fxpxp(6,5,5)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxpxp(7,5,5)=exp(kp);
fxpxp(10,5,5)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)) - PHI*exp(2*g)*exp(-k)*exp(2*kp)*exp(-yy);
fxpxp(6,7,5)=-(BETTA*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxp(6,8,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxp(1,6,6)=(3*PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (2*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)) + (PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fxpxp(4,6,6)=(BETTA*PSSI*exp(2*dp)*exp(-2*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxp(6,5,7)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1));
fxpxp(6,7,7)=(BETTA*(2*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) + (ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxp(6,8,7)=-(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxp(6,5,8)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxp(6,7,8)=-(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpxp(6,8,8)=(ALFA*BETTA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpxp
if any(any(any(isnan(fxpxp),1),2),3) ~= 0 || any(any(any(isinf(fxpxp),1),2),3) ~= 0 || any(any(any(imag(fxpxp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpy
fxpy(6,5,5)=-(GAMA*PHI*THETA*exp(-k)*exp(g)*exp(h)*exp(kp)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxpy(10,5,6)=-PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxpy(1,6,6)=(2*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - (2*PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fxpy(4,6,6)=- (BETTA*PSSI*exp(2*dp)*exp(-2*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpy(6,5,7)=(GAMA*PHI*exp(-k)*exp(c)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fxpy(6,7,9)=-(2*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpy
if any(any(any(isnan(fxpy),1),2),3) ~= 0 || any(any(any(isinf(fxpy),1),2),3) ~= 0 || any(any(any(imag(fxpy),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpyp
fxpyp(6,5,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpyp(4,6,5)=(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyp(6,7,5)=(BETTA*((ALFA^2*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyp(6,8,5)=(ALFA*BETTA*GAMA*THETA*exp(ap)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpyp(6,5,7)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fxpyp(4,6,7)=-(BETTA*GAMA*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyp(6,7,7)=-(BETTA*GAMA*exp(cp)*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyp(6,8,7)=-(ALFA*BETTA*GAMA*exp(ap)*exp(cp)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyp(6,7,9)=(2*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fxpyp
if any(any(any(isnan(fxpyp),1),2),3) ~= 0 || any(any(any(isinf(fxpyp),1),2),3) ~= 0 || any(any(any(imag(fxpyp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyy
fyy(11,1,1)=-exp(gy);
fyy(12,2,2)=-exp(gcc);
fyy(13,3,3)=-exp(giv);
fyy(2,5,5)=(ALFA*exp(2*g)*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^(ALFA + 1) - (exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fyy(4,5,5)=- (GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(5,5,5)=- THETA*exp(h)*exp(h)^(OMEGA - 2)*(OMEGA - 1) - THETA*exp(2*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2) - ALFA*exp(-2*h)*exp(2*k)*exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fyy(6,5,5)=(GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) + (GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) + (GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(4,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(6,7,5)=-(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(1,6,6)=(PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (2*PSSI^2*exp(3*dp)*exp(-2*yy)*exp(g)*exp(2*exp(-yy)*exp(dp) - 2*DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^3 - exp(yy) + (PSSI*exp(3*dp)*exp(-2*yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fyy(2,6,6)=-exp(yy);
fyy(4,6,6)=(BETTA*PSSI*exp(2*dp)*exp(-2*yy)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fyy(10,6,6)=- exp(-yy)*(exp(c) + exp(invs) - exp(yy) + (PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2) - 1;
fyy(11,6,6)=exp(-yyback)*exp(gback)*exp(yy);
fyy(14,6,6)=exp(yy);
fyy(10,7,6)=exp(-yy)*exp(c);
fyy(10,8,6)=exp(-yy)*exp(invs);
fyy(4,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(6,5,7)=-(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(10,6,7)=exp(-yy)*exp(c);
fyy(1,7,7)=exp(c);
fyy(4,7,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*exp(2*c)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(6,7,7)=(GAMA*exp(2*c)*(GAMA + 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*exp(c)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fyy(10,7,7)=-exp(-yy)*exp(c);
fyy(12,7,7)=exp(-cback)*exp(c)*exp(gback);
fyy(15,7,7)=exp(c);
fyy(10,6,8)=exp(-yy)*exp(invs);
fyy(1,8,8)=exp(invs);
fyy(3,8,8)=-exp(invs);
fyy(10,8,8)=-exp(-yy)*exp(invs);
fyy(13,8,8)=exp(-ivback)*exp(gback)*exp(invs);
fyy(16,8,8)=exp(invs);
fyy(6,9,9)=(2*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fyy(7,9,9)=-exp(k1);
% END DISPLAYING fyy
if any(any(any(isnan(fyy),1),2),3) ~= 0 || any(any(any(isinf(fyy),1),2),3) ~= 0 || any(any(any(imag(fyy),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyyp
fyyp(4,6,5)=-(BETTA*GAMA*PSSI*THETA*exp(-yy)*exp(dp)*exp(hp)*exp(exp(-yy)*exp(dp) - DY)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyp(6,9,5)=-(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyp(4,6,7)=(BETTA*GAMA*PSSI*exp(-yy)*exp(cp)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyp(6,9,7)=(BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyyp(6,9,9)=-(2*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fyyp
if any(any(any(isnan(fyyp),1),2),3) ~= 0 || any(any(any(isinf(fyyp),1),2),3) ~= 0 || any(any(any(imag(fyyp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fypyp
fypyp(4,5,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypyp(6,5,5)=(ALFA^2*BETTA*exp(2*gp)*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^(ALFA + 1)) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(gp)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypyp(4,7,5)=-(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypyp(6,7,5)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypyp(6,9,5)=(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypyp(4,5,7)=-(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypyp(6,5,7)=(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fypyp(4,7,7)=(BETTA*GAMA*exp(2*cp)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))*(GAMA + 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*exp(cp)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypyp(6,7,7)=(BETTA*GAMA*exp(cp)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*exp(2*cp)*(GAMA + 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fypyp(6,9,7)=-(BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypyp(6,5,9)=(BETTA*GAMA*PHI*THETA*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(hp)*exp(hp)^(OMEGA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypyp(6,7,9)=-(BETTA*GAMA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p)*exp(cp))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fypyp(6,9,9)=(2*BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fypyp
if any(any(any(isnan(fypyp),1),2),3) ~= 0 || any(any(any(isinf(fypyp),1),2),3) ~= 0 || any(any(any(imag(fypyp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
fxpx  = permute(fxxp,[1,3,2]);
fyx   = permute(fxy,[1,3,2]);
fyxp  = permute(fxpy,[1,3,2]);
fypx  = permute(fxyp,[1,3,2]);
fypxp = permute(fxpyp,[1,3,2]);
fypy  = permute(fyyp,[1,3,2]);
 
if errorMes == 1
    fxx   = NaN;
    fxxp  = NaN;
    fxy   = NaN;
    fxyp  = NaN;
    fxpx  = NaN;
    fxpxp = NaN;
    fxpy  = NaN;
    fxpyp = NaN;
    fyx   = NaN;
    fyxp  = NaN;
    fyy   = NaN;
    fyyp  = NaN;
    fypx  = NaN;
    fypxp = NaN;
    fypy  = NaN;
    fypyp = NaN;
end
end

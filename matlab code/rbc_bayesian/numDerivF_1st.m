function [f,fx,fxp,fy,fyp,eta,errorMes] = numDerivF_1st(params)
 
%Initializing the flag for errors
errorMes = 0;
 
%The steady state of the model
[GAMA, DELTA, ALFA, PSSI, OMEGA, G, SIGMAG, RHOG, SIGMAA, RHOA, PHI, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, la, lap, a, ap, tby, tbyp, g, gp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, gback, gbackp, nx, ny, ne, sig, eta]=growth_model_ss(params);;
 
%Setting dimension for the matrices
n   = ny+nx;
f   = zeros(n,1);
fx  = zeros(n,nx);
fxp = zeros(n,nx);
fy  = zeros(n,ny);
fyp = zeros(n,ny);
 
% START DISPLAYING f
f(1,1)=exp(c) + exp(d) + exp(invs) - exp(yy) + (PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1));
f(2,1)=exp(a)*exp(k)^ALFA*(exp(g)*exp(h))^(1 - ALFA) - exp(yy);
f(3,1)=exp(k)*(DELTA - 1) - exp(invs) + exp(g)*exp(kp);
f(4,1)=(BETTA*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - 1/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
f(5,1)=- THETA*exp(h)^(OMEGA - 1) - exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
f(6,1)=(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1)/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA - (BETTA*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
f(7,1)=exp(kp) - exp(k1);
f(8,1)=RHOA*log(exp(a)) - log(exp(ap));
f(9,1)=RHOG*log(exp(g)/G) - log(exp(gp)/G);
f(10,1)=- tby - exp(-yy)*(exp(c) + exp(invs) - exp(yy) + (PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2);
f(11,1)=exp(-yyback)*exp(gback)*exp(yy) - exp(gy);
f(12,1)=exp(-cback)*exp(c)*exp(gback) - exp(gcc);
f(13,1)=exp(-ivback)*exp(gback)*exp(invs) - exp(giv);
f(14,1)=exp(yy) - exp(yybackp);
f(15,1)=exp(c) - exp(cbackp);
f(16,1)=exp(invs) - exp(ivbackp);
f(17,1)=exp(g) - exp(gbackp);
% END DISPLAYING f
if any(any(isnan(f),1),2) ~= 0 || any(any(isinf(f),1),2) ~= 0 || any(any(imag(f),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fx
fx(11,1)=-exp(-yyback)*exp(gback)*exp(yy);
fx(12,2)=-exp(-cback)*exp(c)*exp(gback);
fx(13,3)=-exp(-ivback)*exp(gback)*exp(invs);
fx(11,4)=exp(-yyback)*exp(gback)*exp(yy);
fx(12,4)=exp(-cback)*exp(c)*exp(gback);
fx(13,4)=exp(-ivback)*exp(gback)*exp(invs);
fx(1,5)=(PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fx(2,5)=ALFA*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(exp(g)*exp(h))^(1 - ALFA);
fx(3,5)=exp(k)*(DELTA - 1);
fx(5,5)=-ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fx(6,5)=(PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fx(10,5)=-exp(-yy)*((PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2 + PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp)));
fx(1,6)=exp(d);
fx(1,7)=- (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)) - PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fx(2,7)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fx(3,7)=exp(g)*exp(kp);
fx(4,7)=-(BETTA*GAMA*exp(g)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fx(5,7)=(exp(a)*exp(g)*(exp(-h)*exp(k))^ALFA*(ALFA - 1)^2)/exp(g)^ALFA;
fx(6,7)=(BETTA*GAMA*exp(g)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^(GAMA + 1)*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fx(9,7)=RHOG;
fx(10,7)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fx(17,7)=exp(g);
fx(2,8)=exp(a)*exp(k)^ALFA*(exp(g)*exp(h))^(1 - ALFA);
fx(5,8)=-exp(a)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
fx(8,8)=RHOA;
% END DISPLAYING fx
if any(any(isnan(fx),1),2) ~= 0 || any(any(isinf(fx),1),2) ~= 0 || any(any(imag(fx),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxp
fxp(14,1)=-exp(yybackp);
fxp(15,2)=-exp(cbackp);
fxp(16,3)=-exp(ivbackp);
fxp(17,4)=-exp(gbackp);
fxp(1,5)=-PHI*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxp(3,5)=exp(g)*exp(kp);
fxp(6,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA) - (PHI*exp(-k)*exp(g)*exp(kp))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
fxp(7,5)=exp(kp);
fxp(10,5)=PHI*exp(-yy)*exp(g)*exp(kp)*(G - exp(-k)*exp(g)*exp(kp));
fxp(1,6)=(PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2 - (exp(dp)*exp(g))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1));
fxp(4,6)=(BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxp(6,7)=(BETTA*(PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p) - (ALFA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(-kp)*exp(gp)*exp(hp))^ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxp(9,7)=-1;
fxp(6,8)=(ALFA*BETTA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxp(8,8)=-1;
% END DISPLAYING fxp
if any(any(isnan(fxp),1),2) ~= 0 || any(any(isinf(fxp),1),2) ~= 0 || any(any(imag(fxp),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fy
fy(11,1)=-exp(gy);
fy(12,2)=-exp(gcc);
fy(13,3)=-exp(giv);
fy(10,4)=-1;
fy(2,5)=-(exp(a)*exp(g)*exp(h)*exp(k)^ALFA*(ALFA - 1))/(exp(g)*exp(h))^ALFA;
fy(4,5)=-(GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(5,5)=ALFA*exp(-h)*exp(a)*exp(k)*exp(g)^(1 - ALFA)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) - THETA*exp(h)*exp(h)^(OMEGA - 2)*(OMEGA - 1);
fy(6,5)=(GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(1,6)=- exp(yy) - (PSSI*exp(2*dp)*exp(-yy)*exp(g)*exp(exp(-yy)*exp(dp) - DY))/(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1))^2;
fy(2,6)=-exp(yy);
fy(4,6)=-(BETTA*PSSI*exp(-yy)*exp(dp)*exp(exp(-yy)*exp(dp) - DY))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fy(10,6)=exp(-yy)*(exp(c) + exp(invs) - exp(yy) + (PHI*exp(k)*(G - exp(-k)*exp(g)*exp(kp))^2)/2) + 1;
fy(11,6)=exp(-yyback)*exp(gback)*exp(yy);
fy(14,6)=exp(yy);
fy(1,7)=exp(c);
fy(4,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(6,7)=-(GAMA*exp(c)*(PHI*(G - exp(-k)*exp(g)*exp(kp)) - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(10,7)=-exp(-yy)*exp(c);
fy(12,7)=exp(-cback)*exp(c)*exp(gback);
fy(15,7)=exp(c);
fy(1,8)=exp(invs);
fy(3,8)=-exp(invs);
fy(10,8)=-exp(-yy)*exp(invs);
fy(13,8)=exp(-ivback)*exp(gback)*exp(invs);
fy(16,8)=exp(invs);
fy(6,9)=-(BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fy(7,9)=-exp(k1);
% END DISPLAYING fy
if any(any(isnan(fy),1),2) ~= 0 || any(any(isinf(fy),1),2) ~= 0 || any(any(imag(fy),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyp
fyp(4,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyp(6,5)=- (BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(gp)*exp(hp)*(ALFA - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(gp)*exp(hp))^ALFA);
fyp(4,7)=-(BETTA*GAMA*exp(cp)*(RSTAR + PSSI*(exp(exp(-yy)*exp(dp) - DY) - 1)))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyp(6,7)=(BETTA*GAMA*exp(cp)*(DELTA + (PHI*(G - exp(-k1)*exp(gp)*exp(k1p))^2)/2 - ALFA*exp(ap)*(exp(-kp)*exp(gp)*exp(hp))^(1 - ALFA) + PHI*exp(-k1)*exp(gp)*exp(k1p)*(G - exp(-k1)*exp(gp)*exp(k1p)) - 1))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fyp(6,9)=(BETTA*PHI*exp(2*gp)*exp(-2*k1)*exp(2*k1p))/(exp(g)^GAMA*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
% END DISPLAYING fyp
if any(any(isnan(fyp),1),2) ~= 0 || any(any(isinf(fyp),1),2) ~= 0 || any(any(imag(fyp),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
if errorMes == 1
    fx  = NaN;
    fxp = NaN;
    fy  = NaN;
    fyp = NaN;
    eta = NaN;
end
end

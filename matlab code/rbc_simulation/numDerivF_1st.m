function [f,fx,fxp,fy,fyp,eta,errorMes] = numDerivF_1st(params)
 
%Initializing the flag for errors
errorMes = 0;
 
%The steady state of the model
[GAMA, DELTA, ALFA, PSSI, OMEGA, SIGMAA, RHOA, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, YY, la, lap, a, ap, tby, tbyp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, nx, ny, ne, sig, eta]=growth_model_ss(params);;
 
%Setting dimension for the matrices
n   = ny+nx;
f   = zeros(n,1);
fx  = zeros(n,nx);
fxp = zeros(n,nx);
fy  = zeros(n,ny);
fyp = zeros(n,ny);
 
% START DISPLAYING f
f(1,1)=exp(c) + exp(d) + exp(invs) - exp(yy) - exp(dp)/(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1));
f(2,1)=exp(a)*exp(h)^(1 - ALFA)*exp(k)^ALFA - exp(yy);
f(3,1)=exp(kp) - exp(invs) + exp(k)*(DELTA - 1);
f(4,1)=(BETTA*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA - 1/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
f(5,1)=- THETA*exp(h)^(OMEGA - 1) - exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
f(6,1)=(BETTA*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA - 1/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^GAMA;
f(7,1)=exp(kp) - exp(k1);
f(8,1)=RHOA*log(exp(a)) - log(exp(ap));
f(9,1)=- tby - exp(-yy)*(exp(c) + exp(invs) - exp(yy));
f(10,1)=exp(-yyback)*exp(yy) - exp(gy);
f(11,1)=exp(-cback)*exp(c) - exp(gcc);
f(12,1)=exp(-ivback)*exp(invs) - exp(giv);
f(13,1)=exp(yy) - exp(yybackp);
f(14,1)=exp(c) - exp(cbackp);
f(15,1)=exp(invs) - exp(ivbackp);
% END DISPLAYING f
if any(any(isnan(f),1),2) ~= 0 || any(any(isinf(f),1),2) ~= 0 || any(any(imag(f),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fx
fx(10,1)=-exp(-yyback)*exp(yy);
fx(11,2)=-exp(-cback)*exp(c);
fx(12,3)=-exp(-ivback)*exp(invs);
fx(2,4)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1);
fx(3,4)=exp(k)*(DELTA - 1);
fx(5,4)=-ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fx(1,5)=exp(d);
fx(2,6)=exp(a)*exp(h)^(1 - ALFA)*exp(k)^ALFA;
fx(5,6)=-exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
fx(8,6)=RHOA;
% END DISPLAYING fx
if any(any(isnan(fx),1),2) ~= 0 || any(any(isinf(fx),1),2) ~= 0 || any(any(imag(fx),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxp
fxp(13,1)=-exp(yybackp);
fxp(14,2)=-exp(cbackp);
fxp(15,3)=-exp(ivbackp);
fxp(3,4)=exp(kp);
fxp(6,4)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxp(7,4)=exp(kp);
fxp(1,5)=(PSSI*exp(2*dp)*exp(exp(dp)/YY - DY))/(YY*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^2) - exp(dp)/(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1));
fxp(4,5)=(BETTA*PSSI*exp(exp(dp)/YY - DY)*exp(dp))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxp(6,6)=(ALFA*BETTA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA;
fxp(8,6)=-1;
% END DISPLAYING fxp
if any(any(isnan(fxp),1),2) ~= 0 || any(any(isinf(fxp),1),2) ~= 0 || any(any(imag(fxp),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fy
fy(10,1)=-exp(gy);
fy(11,2)=-exp(gcc);
fy(12,3)=-exp(giv);
fy(9,4)=-1;
fy(2,5)=-(exp(a)*exp(h)*exp(k)^ALFA*(ALFA - 1))/exp(h)^ALFA;
fy(4,5)=-(GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(5,5)=ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) - THETA*exp(h)*exp(h)^(OMEGA - 2)*(OMEGA - 1);
fy(6,5)=-(GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(1,6)=-exp(yy);
fy(2,6)=-exp(yy);
fy(9,6)=exp(-yy)*(exp(c) + exp(invs) - exp(yy)) + 1;
fy(10,6)=exp(-yyback)*exp(yy);
fy(13,6)=exp(yy);
fy(1,7)=exp(c);
fy(4,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(6,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1);
fy(9,7)=-exp(-yy)*exp(c);
fy(11,7)=exp(-cback)*exp(c);
fy(14,7)=exp(c);
fy(1,8)=exp(invs);
fy(3,8)=-exp(invs);
fy(9,8)=-exp(-yy)*exp(invs);
fy(12,8)=exp(-ivback)*exp(invs);
fy(15,8)=exp(invs);
fy(7,9)=-exp(k1);
% END DISPLAYING fy
if any(any(isnan(fy),1),2) ~= 0 || any(any(isinf(fy),1),2) ~= 0 || any(any(imag(fy),1),2) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyp
fyp(4,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
fyp(6,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fyp(4,7)=-(BETTA*GAMA*exp(cp)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
fyp(6,7)=-(BETTA*GAMA*exp(cp)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
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

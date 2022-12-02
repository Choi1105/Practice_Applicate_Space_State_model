function [fxx,fxxp,fxy,fxyp,fxpx,fxpxp,fxpy,fxpyp,fyx,fyxp,fyy,fyyp,fypx,fypxp,fypy,fypyp,errorMes] = numDerivF_2nd(params)
 
%Initializing the flag for errors
errorMes = 0;
 
%The steady state of the model
[GAMA, DELTA, ALFA, PSSI, OMEGA, SIGMAA, RHOA, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, YY, la, lap, a, ap, tby, tbyp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, nx, ny, ne, sig, eta]=growth_model_ss(params);;
 
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
fxx(10,1,1)=exp(-yyback)*exp(yy);
fxx(11,2,2)=exp(-cback)*exp(c);
fxx(12,3,3)=exp(-ivback)*exp(invs);
fxx(2,4,4)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1) + ALFA*exp(2*k)*exp(a)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 2)*(ALFA - 1);
fxx(3,4,4)=exp(k)*(DELTA - 1);
fxx(5,4,4)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxx(2,6,4)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1);
fxx(5,6,4)=-ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxx(1,5,5)=exp(d);
fxx(2,4,6)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1);
fxx(5,4,6)=-ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxx(2,6,6)=exp(a)*exp(h)^(1 - ALFA)*exp(k)^ALFA;
fxx(5,6,6)=-exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
% END DISPLAYING fxx
if any(any(any(isnan(fxx),1),2),3) ~= 0 || any(any(any(isinf(fxx),1),2),3) ~= 0 || any(any(any(imag(fxx),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxp
% END DISPLAYING fxxp
if any(any(any(isnan(fxxp),1),2),3) ~= 0 || any(any(any(isinf(fxxp),1),2),3) ~= 0 || any(any(any(imag(fxxp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxy
fxy(2,4,5)=-(ALFA*exp(a)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/exp(h)^ALFA;
fxy(5,4,5)=ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxy(2,6,5)=-(exp(a)*exp(h)*exp(k)^ALFA*(ALFA - 1))/exp(h)^ALFA;
fxy(5,6,5)=ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxy(10,1,6)=-exp(-yyback)*exp(yy);
fxy(11,2,7)=-exp(-cback)*exp(c);
fxy(12,3,8)=-exp(-ivback)*exp(invs);
% END DISPLAYING fxy
if any(any(any(isnan(fxy),1),2),3) ~= 0 || any(any(any(isinf(fxy),1),2),3) ~= 0 || any(any(any(imag(fxy),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxyp
% END DISPLAYING fxyp
if any(any(any(isnan(fxyp),1),2),3) ~= 0 || any(any(any(isinf(fxyp),1),2),3) ~= 0 || any(any(any(imag(fxyp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxp
fxpxp(13,1,1)=-exp(yybackp);
fxpxp(14,2,2)=-exp(cbackp);
fxpxp(15,3,3)=-exp(ivbackp);
fxpxp(3,4,4)=exp(kp);
fxpxp(6,4,4)=(ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxp(7,4,4)=exp(kp);
fxpxp(6,6,4)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxp(1,5,5)=(3*PSSI*exp(2*dp)*exp(exp(dp)/YY - DY))/(YY*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^2) - exp(dp)/(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)) + (PSSI*exp(3*dp)*exp(exp(dp)/YY - DY))/(YY^2*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^2) - (2*PSSI^2*exp(3*dp)*exp((2*exp(dp))/YY - 2*DY))/(YY^2*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^3);
fxpxp(4,5,5)=(BETTA*PSSI*exp(2*dp)*exp(exp(dp)/YY - DY))/(YY^2*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(exp(dp)/YY - DY)*exp(dp))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxp(6,4,6)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxp(6,6,6)=(ALFA*BETTA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA;
% END DISPLAYING fxpxp
if any(any(any(isnan(fxpxp),1),2),3) ~= 0 || any(any(any(isinf(fxpxp),1),2),3) ~= 0 || any(any(any(imag(fxpxp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpy
% END DISPLAYING fxpy
if any(any(any(isnan(fxpy),1),2),3) ~= 0 || any(any(any(isinf(fxpy),1),2),3) ~= 0 || any(any(any(imag(fxpy),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpyp
fxpyp(6,4,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpyp(4,5,5)=(BETTA*GAMA*PSSI*THETA*exp(exp(dp)/YY - DY)*exp(dp)*exp(hp)*exp(hp)^(OMEGA - 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyp(6,6,5)=(ALFA*BETTA*GAMA*THETA*exp(ap)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpyp(6,4,7)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpyp(4,5,7)=-(BETTA*GAMA*PSSI*exp(exp(dp)/YY - DY)*exp(cp)*exp(dp))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpyp(6,6,7)=-(ALFA*BETTA*GAMA*exp(ap)*exp(cp)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
% END DISPLAYING fxpyp
if any(any(any(isnan(fxpyp),1),2),3) ~= 0 || any(any(any(isinf(fxpyp),1),2),3) ~= 0 || any(any(any(imag(fxpyp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyy
fyy(10,1,1)=-exp(gy);
fyy(11,2,2)=-exp(gcc);
fyy(12,3,3)=-exp(giv);
fyy(2,5,5)=(ALFA*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/exp(h)^(ALFA + 1) - (exp(a)*exp(h)*exp(k)^ALFA*(ALFA - 1))/exp(h)^ALFA;
fyy(4,5,5)=- (GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(5,5,5)=- THETA*exp(h)*exp(h)^(OMEGA - 2)*(OMEGA - 1) - THETA*exp(2*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2) - ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fyy(6,5,5)=- (GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(4,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(6,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(1,6,6)=-exp(yy);
fyy(2,6,6)=-exp(yy);
fyy(9,6,6)=- exp(-yy)*(exp(c) + exp(invs) - exp(yy)) - 1;
fyy(10,6,6)=exp(-yyback)*exp(yy);
fyy(13,6,6)=exp(yy);
fyy(9,7,6)=exp(-yy)*exp(c);
fyy(9,8,6)=exp(-yy)*exp(invs);
fyy(4,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(6,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(9,6,7)=exp(-yy)*exp(c);
fyy(1,7,7)=exp(c);
fyy(4,7,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*exp(2*c)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(6,7,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*exp(2*c)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyy(9,7,7)=-exp(-yy)*exp(c);
fyy(11,7,7)=exp(-cback)*exp(c);
fyy(14,7,7)=exp(c);
fyy(9,6,8)=exp(-yy)*exp(invs);
fyy(1,8,8)=exp(invs);
fyy(3,8,8)=-exp(invs);
fyy(9,8,8)=-exp(-yy)*exp(invs);
fyy(12,8,8)=exp(-ivback)*exp(invs);
fyy(15,8,8)=exp(invs);
fyy(7,9,9)=-exp(k1);
% END DISPLAYING fyy
if any(any(any(isnan(fyy),1),2),3) ~= 0 || any(any(any(isinf(fyy),1),2),3) ~= 0 || any(any(any(imag(fyy),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyyp
% END DISPLAYING fyyp
if any(any(any(isnan(fyyp),1),2),3) ~= 0 || any(any(any(isinf(fyyp),1),2),3) ~= 0 || any(any(any(imag(fyyp),1),2),3) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fypyp
fypyp(4,5,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(OMEGA - 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
fypyp(6,5,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) + (ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) - (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fypyp(4,7,5)=-(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypyp(6,7,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypyp(4,5,7)=-(BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypyp(6,5,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypyp(4,7,7)=(BETTA*GAMA*exp(2*cp)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (BETTA*GAMA*exp(cp)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
fypyp(6,7,7)=(BETTA*GAMA*exp(2*cp)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (BETTA*GAMA*exp(cp)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
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

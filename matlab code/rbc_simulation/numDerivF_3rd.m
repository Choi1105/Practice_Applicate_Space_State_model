function [fxxx,fxxxp,fxxy,fxxyp,fxxpxp,fxxpy,fxxpyp,fxyy,fxyyp,fxypyp,fxpxpxp,fxpxpy,fxpxpyp,fxpyy,fxpyyp,fxpypyp,fyyy,fyyyp,fyypyp,fypypyp,errorMes] = numDerivF_3rd(params)
%Initializing the flag for errors
errorMes = 0;
 
%The steady state of the model
[GAMA, DELTA, ALFA, PSSI, OMEGA, SIGMAA, RHOA, RSTAR, BETTA, THETA, c, cp, h, hp, k, kp, k1, k1p, d, dp, invs, invsp, tb, tbp, DY, YY, la, lap, a, ap, tby, tbyp, gy, gyp, gcc, gccp, giv, givp, yy, yyp, yyback, yybackp, cback, cbackp, ivback, ivbackp, nx, ny, ne, sig, eta]=growth_model_ss(params);;
 
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
fxxx(10,1,1,1)=-exp(-yyback)*exp(yy);
fxxx(11,2,2,2)=-exp(-cback)*exp(c);
fxxx(12,3,3,3)=-exp(-ivback)*exp(invs);
fxxx(2,4,4,4)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1) + 3*ALFA*exp(2*k)*exp(a)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 2)*(ALFA - 1) + ALFA*exp(3*k)*exp(a)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 3)*(ALFA - 1)*(ALFA - 2);
fxxx(3,4,4,4)=exp(k)*(DELTA - 1);
fxxx(5,4,4,4)=- 3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) - ALFA*exp(-3*h)*exp(3*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fxxx(2,6,4,4)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1) + ALFA*exp(2*k)*exp(a)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 2)*(ALFA - 1);
fxxx(5,6,4,4)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,4,6,4)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1) + ALFA*exp(2*k)*exp(a)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 2)*(ALFA - 1);
fxxx(5,4,6,4)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,6,6,4)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1);
fxxx(5,6,6,4)=-ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(1,5,5,5)=exp(d);
fxxx(2,4,4,6)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1) + ALFA*exp(2*k)*exp(a)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 2)*(ALFA - 1);
fxxx(5,4,4,6)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,6,4,6)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1);
fxxx(5,6,4,6)=-ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,4,6,6)=ALFA*exp(a)*exp(k)*exp(h)^(1 - ALFA)*exp(k)^(ALFA - 1);
fxxx(5,4,6,6)=-ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxx(2,6,6,6)=exp(a)*exp(h)^(1 - ALFA)*exp(k)^ALFA;
fxxx(5,6,6,6)=-exp(a)*(exp(-h)*exp(k))^ALFA*(ALFA - 1);
% END DISPLAYING fxxx
if any(any(any(any(isnan(fxxx),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxx),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxx),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxxp
% END DISPLAYING fxxxp
if any(any(any(any(isnan(fxxxp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxxp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxxp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxy
fxxy(2,4,4,5)=- (ALFA*exp(2*k)*exp(a)*exp(h)*exp(k)^(ALFA - 2)*(ALFA - 1)^2)/exp(h)^ALFA - (ALFA*exp(a)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/exp(h)^ALFA;
fxxy(5,4,4,5)=3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) + ALFA*exp(-3*h)*exp(3*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fxxy(2,6,4,5)=-(ALFA*exp(a)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/exp(h)^ALFA;
fxxy(5,6,4,5)=ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxy(2,4,6,5)=-(ALFA*exp(a)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/exp(h)^ALFA;
fxxy(5,4,6,5)=ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 + ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxy(2,6,6,5)=-(exp(a)*exp(h)*exp(k)^ALFA*(ALFA - 1))/exp(h)^ALFA;
fxxy(5,6,6,5)=ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxxy(10,1,1,6)=exp(-yyback)*exp(yy);
fxxy(11,2,2,7)=exp(-cback)*exp(c);
fxxy(12,3,3,8)=exp(-ivback)*exp(invs);
% END DISPLAYING fxxy
if any(any(any(any(isnan(fxxy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxyp
% END DISPLAYING fxxyp
if any(any(any(any(isnan(fxxyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxpxp
% END DISPLAYING fxxpxp
if any(any(any(any(isnan(fxxpxp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxpxp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxpxp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxpy
% END DISPLAYING fxxpy
if any(any(any(any(isnan(fxxpy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxpy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxpy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxxpyp
% END DISPLAYING fxxpyp
if any(any(any(any(isnan(fxxpyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxxpyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxxpyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxyy
fxyy(2,4,5,5)=(ALFA^2*exp(2*h)*exp(a)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/exp(h)^(ALFA + 1) - (ALFA*exp(a)*exp(h)*exp(k)*exp(k)^(ALFA - 1)*(ALFA - 1))/exp(h)^ALFA;
fxyy(5,4,5,5)=- 3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) - ALFA*exp(-3*h)*exp(3*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fxyy(2,6,5,5)=(ALFA*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/exp(h)^(ALFA + 1) - (exp(a)*exp(h)*exp(k)^ALFA*(ALFA - 1))/exp(h)^ALFA;
fxyy(5,6,5,5)=- ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1);
fxyy(10,1,6,6)=-exp(-yyback)*exp(yy);
fxyy(11,2,7,7)=-exp(-cback)*exp(c);
fxyy(12,3,8,8)=-exp(-ivback)*exp(invs);
% END DISPLAYING fxyy
if any(any(any(any(isnan(fxyy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxyy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxyy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxyyp
% END DISPLAYING fxyyp
if any(any(any(any(isnan(fxyyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxyyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxyyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxypyp
% END DISPLAYING fxypyp
if any(any(any(any(isnan(fxypyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxypyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxypyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxpxp
fxpxpxp(13,1,1,1)=-exp(yybackp);
fxpxpxp(14,2,2,2)=-exp(cbackp);
fxpxpxp(15,3,3,3)=-exp(ivbackp);
fxpxpxp(3,4,4,4)=exp(kp);
fxpxpxp(6,4,4,4)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) - (3*ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) + (ALFA^2*BETTA*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 2));
fxpxpxp(7,4,4,4)=exp(kp);
fxpxpxp(6,6,4,4)=(ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxpxp(6,4,6,4)=(ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxpxp(6,6,6,4)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxpxp(1,5,5,5)=(7*PSSI*exp(2*dp)*exp(exp(dp)/YY - DY))/(YY*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^2) - exp(dp)/(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)) + (6*PSSI*exp(3*dp)*exp(exp(dp)/YY - DY))/(YY^2*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^2) + (PSSI*exp(4*dp)*exp(exp(dp)/YY - DY))/(YY^3*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^2) - (12*PSSI^2*exp(3*dp)*exp((2*exp(dp))/YY - 2*DY))/(YY^2*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^3) - (6*PSSI^2*exp(4*dp)*exp((2*exp(dp))/YY - 2*DY))/(YY^3*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^3) + (6*PSSI^3*exp(4*dp)*exp(exp(dp)/YY - DY)*exp((2*exp(dp))/YY - 2*DY))/(YY^3*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))^4);
fxpxpxp(4,5,5,5)=(3*BETTA*PSSI*exp(2*dp)*exp(exp(dp)/YY - DY))/(YY^2*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(3*dp)*exp(exp(dp)/YY - DY))/(YY^3*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA) + (BETTA*PSSI*exp(exp(dp)/YY - DY)*exp(dp))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA);
fxpxpxp(6,4,4,6)=(ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxpxp(6,6,4,6)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxpxp(6,4,6,6)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxpxp(6,6,6,6)=(ALFA*BETTA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA;
% END DISPLAYING fxpxpxp
if any(any(any(any(isnan(fxpxpxp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpxpxp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpxpxp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxpy
% END DISPLAYING fxpxpy
if any(any(any(any(isnan(fxpxpy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpxpy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpxpy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpxpyp
fxpxpyp(6,4,4,5)=(3*ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 2)) + (ALFA^2*BETTA*GAMA*THETA*exp(3*hp)*exp(-2*kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpxpyp(6,6,4,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpxpyp(4,5,5,5)=(BETTA*GAMA*PSSI*THETA*exp(exp(dp)/YY - DY)*exp(dp)*exp(hp)*exp(hp)^(OMEGA - 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PSSI*THETA*exp(2*dp)*exp(exp(dp)/YY - DY)*exp(hp)*exp(hp)^(OMEGA - 1))/(YY^2*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,4,6,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) - (ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpxpyp(6,6,6,5)=(ALFA*BETTA*GAMA*THETA*exp(ap)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA);
fxpxpyp(6,4,4,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (ALFA^2*BETTA*GAMA*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1));
fxpxpyp(6,6,4,7)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpxpyp(4,5,5,7)=- (BETTA*GAMA*PSSI*exp(exp(dp)/YY - DY)*exp(cp)*exp(dp))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) - (BETTA*GAMA*PSSI*exp(2*dp)*exp(exp(dp)/YY - DY)*exp(cp))/(YY^2*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpxpyp(6,4,6,7)=-(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpxpyp(6,6,6,7)=-(ALFA*BETTA*GAMA*exp(ap)*exp(cp)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
% END DISPLAYING fxpxpyp
if any(any(any(any(isnan(fxpxpyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpxpyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpxpyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpyy
% END DISPLAYING fxpyy
if any(any(any(any(isnan(fxpyy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpyy),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpyy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpyyp
% END DISPLAYING fxpyyp
if any(any(any(any(isnan(fxpyyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpyyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpyyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fxpypyp
fxpypyp(6,4,5,5)=(ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) - (3*ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) + (ALFA^2*BETTA*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 2)) - (2*ALFA^2*BETTA*GAMA*THETA*exp(3*hp)*exp(-2*kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) + (3*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) + (ALFA*BETTA*GAMA*THETA^2*exp(3*hp)*exp(-kp)*exp(ap)*exp(hp)^(2*OMEGA - 2)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA) + (ALFA*BETTA*GAMA*THETA*exp(3*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 2)*(ALFA - 1)*(OMEGA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpypyp(4,5,5,5)=(BETTA*GAMA*PSSI*THETA*exp(exp(dp)/YY - DY)*exp(dp)*exp(hp)*exp(hp)^(OMEGA - 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)) + (BETTA*GAMA*PSSI*THETA^2*exp(2*hp)*exp(exp(dp)/YY - DY)*exp(dp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) + (BETTA*GAMA*PSSI*THETA*exp(2*hp)*exp(exp(dp)/YY - DY)*exp(dp)*exp(hp)^(OMEGA - 2)*(OMEGA - 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,6,5,5)=(ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) + (ALFA*BETTA*GAMA*THETA*exp(ap)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(ap)*exp(hp)^(OMEGA - 2)*(exp(-kp)*exp(hp))^(1 - ALFA)*(OMEGA - 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (ALFA*BETTA*GAMA*THETA^2*exp(2*hp)*exp(ap)*exp(hp)^(2*OMEGA - 2)*(exp(-kp)*exp(hp))^(1 - ALFA)*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpypyp(6,4,7,5)=(ALFA^2*BETTA*GAMA*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fxpypyp(4,5,7,5)=-(BETTA*GAMA*PSSI*THETA*exp(exp(dp)/YY - DY)*exp(cp)*exp(dp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,6,7,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(ap)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(hp))^(1 - ALFA)*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fxpypyp(6,4,5,7)=(ALFA^2*BETTA*GAMA*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fxpypyp(4,5,5,7)=-(BETTA*GAMA*PSSI*THETA*exp(exp(dp)/YY - DY)*exp(cp)*exp(dp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2));
fxpypyp(6,6,5,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*THETA*exp(ap)*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(exp(-kp)*exp(hp))^(1 - ALFA)*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fxpypyp(6,4,7,7)=(ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA) - (ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fxpypyp(4,5,7,7)=(BETTA*GAMA*PSSI*exp(2*cp)*exp(exp(dp)/YY - DY)*exp(dp)*(GAMA + 1))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)) - (BETTA*GAMA*PSSI*exp(exp(dp)/YY - DY)*exp(cp)*exp(dp))/(YY*(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1));
fxpypyp(6,6,7,7)=(ALFA*BETTA*GAMA*exp(2*cp)*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA)*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (ALFA*BETTA*GAMA*exp(ap)*exp(cp)*(exp(-kp)*exp(hp))^(1 - ALFA))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1);
% END DISPLAYING fxpypyp
if any(any(any(any(isnan(fxpypyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fxpypyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fxpypyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyyy
fyyy(10,1,1,1)=-exp(gy);
fyyy(11,2,2,2)=-exp(gcc);
fyyy(12,3,3,3)=-exp(giv);
fyyy(2,5,5,5)=(3*ALFA*exp(2*h)*exp(a)*exp(k)^ALFA*(ALFA - 1))/exp(h)^(ALFA + 1) - (exp(a)*exp(h)*exp(k)^ALFA*(ALFA - 1))/exp(h)^ALFA - (ALFA*exp(3*h)*exp(a)*exp(k)^ALFA*(ALFA - 1)*(ALFA + 1))/exp(h)^(ALFA + 2);
fyyy(4,5,5,5)=- (GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(3*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA^2*exp(3*h)*exp(h)^(2*OMEGA - 3)*(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^2*exp(3*h)*exp(h)^(OMEGA - 1)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^3*exp(3*h)*exp(h)^(2*OMEGA - 2)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(5,5,5,5)=3*ALFA*exp(-2*h)*exp(2*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 2)*(ALFA - 1)^2 - 3*THETA*exp(2*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2) - THETA*exp(h)*exp(h)^(OMEGA - 2)*(OMEGA - 1) + ALFA*exp(-h)*exp(a)*exp(k)*(exp(-h)*exp(k))^(ALFA - 1)*(ALFA - 1) - THETA*exp(3*h)*exp(h)^(OMEGA - 4)*(OMEGA - 1)*(OMEGA - 2)*(OMEGA - 3) + ALFA*exp(-3*h)*exp(3*k)*exp(a)*(exp(-h)*exp(k))^(ALFA - 3)*(ALFA - 1)^2*(ALFA - 2);
fyyy(6,5,5,5)=- (GAMA*THETA*exp(h)*exp(h)^(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*THETA*exp(2*h)*exp(h)^(OMEGA - 2)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*THETA^2*exp(2*h)*exp(h)^(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(3*h)*exp(h)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (GAMA*THETA^2*exp(3*h)*exp(h)^(2*OMEGA - 3)*(2*OMEGA - 2)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^2*exp(3*h)*exp(h)^(OMEGA - 1)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA^3*exp(3*h)*exp(h)^(2*OMEGA - 2)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(4,7,5,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(6,7,5,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(4,5,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(6,5,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(4,7,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,7,7,5)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(1,6,6,6)=-exp(yy);
fyyy(2,6,6,6)=-exp(yy);
fyyy(9,6,6,6)=exp(-yy)*(exp(c) + exp(invs) - exp(yy)) + 1;
fyyy(10,6,6,6)=exp(-yyback)*exp(yy);
fyyy(13,6,6,6)=exp(yy);
fyyy(9,7,6,6)=-exp(-yy)*exp(c);
fyyy(9,8,6,6)=-exp(-yy)*exp(invs);
fyyy(9,6,7,6)=-exp(-yy)*exp(c);
fyyy(9,7,7,6)=exp(-yy)*exp(c);
fyyy(9,6,8,6)=-exp(-yy)*exp(invs);
fyyy(9,8,8,6)=exp(-yy)*exp(invs);
fyyy(4,5,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(6,5,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*THETA^2*exp(2*h)*exp(c)*exp(h)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3) + (GAMA*THETA*exp(2*h)*exp(c)*exp(h)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2);
fyyy(4,7,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,7,5,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(9,6,6,7)=-exp(-yy)*exp(c);
fyyy(9,7,6,7)=exp(-yy)*exp(c);
fyyy(4,5,7,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,5,7,7)=(GAMA*THETA*exp(c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) - (GAMA*THETA*exp(2*c)*exp(h)*exp(h)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(9,6,7,7)=exp(-yy)*exp(c);
fyyy(1,7,7,7)=exp(c);
fyyy(4,7,7,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*exp(2*c)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*exp(3*c)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(6,7,7,7)=(GAMA*exp(c))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 1) - (3*GAMA*exp(2*c)*(GAMA + 1))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 2) + (GAMA*exp(3*c)*(GAMA + 1)*(GAMA + 2))/(exp(c) - (THETA*exp(h)^OMEGA)/OMEGA)^(GAMA + 3);
fyyy(9,7,7,7)=-exp(-yy)*exp(c);
fyyy(11,7,7,7)=exp(-cback)*exp(c);
fyyy(14,7,7,7)=exp(c);
fyyy(9,6,6,8)=-exp(-yy)*exp(invs);
fyyy(9,8,6,8)=exp(-yy)*exp(invs);
fyyy(9,6,8,8)=exp(-yy)*exp(invs);
fyyy(1,8,8,8)=exp(invs);
fyyy(3,8,8,8)=-exp(invs);
fyyy(9,8,8,8)=-exp(-yy)*exp(invs);
fyyy(12,8,8,8)=exp(-ivback)*exp(invs);
fyyy(15,8,8,8)=exp(invs);
fyyy(7,9,9,9)=-exp(k1);
% END DISPLAYING fyyy
if any(any(any(any(isnan(fyyy),1),2),3),4) ~= 0 || any(any(any(any(isinf(fyyy),1),2),3),4) ~= 0 || any(any(any(any(imag(fyyy),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyyyp
% END DISPLAYING fyyyp
if any(any(any(any(isnan(fyyyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fyyyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fyyyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fyypyp
% END DISPLAYING fyypyp
if any(any(any(any(isnan(fyypyp),1),2),3),4) ~= 0 || any(any(any(any(isinf(fyypyp),1),2),3),4) ~= 0 || any(any(any(any(imag(fyypyp),1),2),3),4) ~= 0
   errorMes = max(errorMes,1);
end
 
% START DISPLAYING fypypyp
fypypyp(4,5,5,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (3*BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (3*BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(OMEGA - 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (BETTA*GAMA*THETA*exp(3*hp)*exp(hp)^(OMEGA - 3)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(OMEGA - 1)*(OMEGA - 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(2*OMEGA - 3)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(2*OMEGA - 2)*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(OMEGA - 1)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA^3*exp(3*hp)*exp(hp)^(2*OMEGA - 2)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3);
fypypyp(6,5,5,5)=(BETTA*GAMA*THETA*exp(hp)*exp(hp)^(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (3*BETTA*GAMA*THETA^2*exp(2*hp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (3*BETTA*GAMA*THETA*exp(2*hp)*exp(hp)^(OMEGA - 2)*(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) - (ALFA*BETTA*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^ALFA) + (3*ALFA^2*BETTA*exp(2*hp)*exp(-2*kp)*exp(ap)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 1)) + (BETTA*GAMA*THETA*exp(3*hp)*exp(hp)^(OMEGA - 3)*(OMEGA - 1)*(OMEGA - 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) + (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(2*OMEGA - 3)*(2*OMEGA - 2)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (ALFA^2*BETTA*exp(3*hp)*exp(-3*kp)*exp(ap)*(ALFA - 1)*(ALFA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^GAMA*(exp(-kp)*exp(hp))^(ALFA + 2)) + (BETTA*GAMA*THETA^2*exp(3*hp)*exp(hp)^(OMEGA - 1)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA^3*exp(3*hp)*exp(hp)^(2*OMEGA - 2)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) + (3*ALFA^2*BETTA*GAMA*THETA*exp(3*hp)*exp(-2*kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) - (6*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 1)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (3*ALFA*BETTA*GAMA*THETA^2*exp(3*hp)*exp(-kp)*exp(ap)*exp(hp)^(2*OMEGA - 2)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA) - (3*ALFA*BETTA*GAMA*THETA*exp(3*hp)*exp(-kp)*exp(ap)*exp(hp)^(OMEGA - 2)*(ALFA - 1)*(OMEGA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA);
fypypyp(4,7,5,5)=- (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypypyp(6,7,5,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (ALFA^2*BETTA*GAMA*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fypypyp(4,5,7,5)=- (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypypyp(6,5,7,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (ALFA^2*BETTA*GAMA*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fypypyp(4,7,7,5)=(BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypypyp(6,7,7,5)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fypypyp(4,5,5,7)=- (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(OMEGA - 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypypyp(6,5,5,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (ALFA^2*BETTA*GAMA*exp(2*hp)*exp(-2*kp)*exp(ap)*exp(cp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^(ALFA + 1)) - (BETTA*GAMA*THETA^2*exp(2*hp)*exp(cp)*exp(hp)^(2*OMEGA - 2)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(2*hp)*exp(cp)*exp(hp)^(OMEGA - 2)*(GAMA + 1)*(OMEGA - 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (2*ALFA*BETTA*GAMA*THETA*exp(2*hp)*exp(-kp)*exp(ap)*exp(cp)*exp(hp)^(OMEGA - 1)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fypypyp(4,7,5,7)=(BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypypyp(6,7,5,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fypypyp(4,5,7,7)=(BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2);
fypypyp(6,5,7,7)=(ALFA*BETTA*GAMA*exp(-kp)*exp(ap)*exp(cp)*exp(hp)*(ALFA - 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1)*(exp(-kp)*exp(hp))^ALFA) - (BETTA*GAMA*THETA*exp(cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) + (BETTA*GAMA*THETA*exp(2*cp)*exp(hp)*exp(hp)^(OMEGA - 1)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3) - (ALFA*BETTA*GAMA*exp(2*cp)*exp(-kp)*exp(ap)*exp(hp)*(ALFA - 1)*(GAMA + 1))/((exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2)*(exp(-kp)*exp(hp))^ALFA);
fypypyp(4,7,7,7)=(3*BETTA*GAMA*exp(2*cp)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (BETTA*GAMA*exp(cp)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1)))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) - (BETTA*GAMA*exp(3*cp)*(RSTAR + PSSI*(exp(exp(dp)/YY - DY) - 1))*(GAMA + 1)*(GAMA + 2))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3);
fypypyp(6,7,7,7)=(3*BETTA*GAMA*exp(2*cp)*(GAMA + 1)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 2) - (BETTA*GAMA*exp(cp)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 1) - (BETTA*GAMA*exp(3*cp)*(GAMA + 1)*(GAMA + 2)*(ALFA*exp(ap)*(exp(-kp)*exp(hp))^(1 - ALFA) - DELTA + 1))/(exp(cp) - (THETA*exp(hp)^OMEGA)/OMEGA)^(GAMA + 3);
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

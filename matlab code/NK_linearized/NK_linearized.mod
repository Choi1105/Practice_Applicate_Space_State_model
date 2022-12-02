var pi x rf i yf y a g r; 

varexo ea eg er;

parameters sigma chi beta phi zeta psi rhoa rhog rhor phipi phiy gamma;

load param_dnk_reduce;
set_param_value('beta',beta);
set_param_value('sigma',sigma);
set_param_value('chi',chi);
set_param_value('psi',psi);
set_param_value('phi',phi);
set_param_value('zeta',zeta);
set_param_value('rhoa',rhoa);
set_param_value('rhog',rhog);
set_param_value('rhor',rhor);
set_param_value('phipi',phipi);
set_param_value('phiy',phiy);
set_param_value('gamma',gamma);

model(linear);

% (1) Phillips Curve
pi = zeta*gamma*x + beta*pi(+1);

% (23) IS equation
x = x(+1) - ((1-psi)/sigma)*(i - pi(+1) - rf);

% (3) rf
rf = (sigma/(1-psi)) * ((1+chi)*(1-psi)/(chi*(1-psi) + sigma))*(a(+1) - a) - (sigma/(1-psi)) * (psi*chi*(1-psi)/(chi*(1-psi)+sigma))*(g(+1) - g);

% (4) yf
yf = ( (1+chi)*(1-psi) / (chi*(1-psi)+sigma))*a + (psi*sigma/(chi*(1-psi) + sigma))*g;

% (5) Output gap
x = y - yf;

% (6) Taylor rule
i = rhor*i(-1) + (1-rhor)*(phipi*pi +  phiy*x) + er;

% (7) Productivity process
a = rhoa*a(-1) + ea;

% (8) Government spending process
g = rhog*g(-1) + eg;

% (9) Fisher relationship
r = i - pi(+1);

end;

shocks;
var ea = 1;
var eg = 1;
var er = 1;
end;

stoch_simul(order=1,irf=20,nograph,ar=1);
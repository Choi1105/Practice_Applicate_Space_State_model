var 
y 
i 
y_f
real_i
n 
a 
pi ;

varexo ea ei;

parameters beta phi phi_pi chi i_star pi_star rho_a rho_i; 

beta = 0.99;
%phi = 0.75;
%phi = 0.9;
phi = 0.5;
phi_pi = 1.5;
chi = 1;
i_star = 0.025;
pi_star = 0;
rho_a = 0.95;
rho_i = 0.8;

model(linear);

% Output function:
y = y(+1)-(i-pi(+1));

% Phillips curve:
pi = ((1-phi)/phi)*(1-phi*beta)*((1+chi)*y-y_f)+beta*pi(+1);

% Flexible output:
y_f = a;

% Productivity shock:
a = rho_a*a(-1) + ea;

% Taylor rule(ignore output gap):
i = (1-rho_i)*i_star+rho_i*i(-1) + (1-rho_i)*phi_pi*(pi-pi_star) + ei;

% Real interest
i = real_i + pi;

% Employment
n = y - a;


end;


shocks;
var ea = 0.01;
var ei = 0.0025;

end;


stoch_simul(order=1,irf=20) y n i real_i pi;

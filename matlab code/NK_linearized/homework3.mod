var y i pi a;

varexo ea ei;

parameters sigma phi beta eta theta_pi rho_a rho_i; 

sigma = 1;
phi = 0.75;
beta = 0.99;
eta = 0.5;
theta_pi = 1.5;
rho_a = 0.9;
rho_i = 0.7;

model(linear);

% Output function:
y = y(+1)-(1/sigma)*(i-pi(+1));

% Phillips curve:
pi = ((1-phi)/phi)*(1-phi*beta)*((eta+sigma)*y-(1+eta)*a)+beta*pi(+1);

% Productivity shock:
a = rho_a*a(-1) + ea;

% Taylor rule(ignore output gap):
i = (1-rho_i)*(1/beta-1) + rho_i*i(-1) + (1-rho_i)*theta_pi*pi + ei;

end;


shocks;
var ea = 1;
var ei = 1;

end;


stoch_simul(order=1,irf=20);



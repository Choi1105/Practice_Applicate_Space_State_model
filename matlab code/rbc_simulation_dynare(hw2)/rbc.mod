// RBC model

var c
    k
    a
    n
    y
    u
    invest
    mu_c
; 

predetermined_variables k;

%Define parameters
parameters delta
        alpha
        phi
        rho_a
        rho_u
        sig_a
        sig_u
        beta
        psi
;

varexo eps_a eps_u
;

rho_a = 0.95;
rho_u = 0.95;
sig_a = 0.01;
sig_u = 0.01;
delta = 0.025; %Depreciation rate
alpha = 0.33; %Capital elasticity of the production function
phi   = 1;
beta  = 0.99;
psi   = 2;
        
model;
%1. Marginal utility of consumption
mu_c = 1/c;

%2 FOC of Labor
mu_c = (psi*n^(phi+alpha))/((1-alpha)*a*k^alpha);

%3 FOC of capital
mu_c = u*beta*mu_c(+1)*(alpha*a(+1)*(k(+1)^(alpha-1)))*((n(+1))^(1-alpha)) + (1/u(+1)*(1-delta));

%4 Resource constraint (see the remark on the fixed typo in the preamble)
y= c + invest;

%5 Definition output
y= a*k^alpha*(n)^(1-alpha);

%6 Definition investment
invest = 1/u*(k(+1) - (1-delta)*k);

%7 LOM temporary TFP
log(a)=rho_a * log(a(-1)) +eps_a;

%8 LOM Invest TFP 
log(u)=rho_u * log(u(-1)) +eps_u;

end;

steady_state_model;
    a = 1;
    u = 1;
    k_n = (((1/beta)-(1-delta))/alpha)^(1/(alpha-1));
    n = 0.3;
    k = k_n*n;
    invest = (delta/u)*k;
    y = (k_n^alpha)*n;
    c   = (k/n)^alpha*n-delta*k;
    mu_c = 1/c;
end;

steady;
check;

shocks;
    var eps_a; stderr 0.01;
    var eps_u; stderr 0.01;
end;


stoch_simul(order=1,irf=20);
// RBC model

var c
    k
    a
    h
    d
    y
    invest
    tb
    mu_c
    tb_y
    g_y
    g_c
    g_invest
    r
; 

predetermined_variables k d;

%Define parameters
parameters gamma
        delta
        alpha
        psi
        omega
        theta
        rho_a
        RSTAR
        beta
;

varexo eps_a
;

rho_a = 0.1;


gamma = 2; %intertemporal elasticity of substitution
delta = 1.03^4-1;%0.03; %Depreciation rate
alpha = 0.32; %Capital elasticity of the production function
omega = 1.6; %exponent of labor in utility function
theta = 1.4*omega;
psi   = 0.001;
RSTAR = 1.1;
beta  = 1/RSTAR;
        
model;
#ru   = RSTAR; %World interest rate
#tb_yu=0.003;
#k_over_ghu =((1/beta-1+delta)/alpha)^(1/(alpha-1)); %k/(g*h)
#hu   = ((1-alpha)*k_over_ghu^alpha/theta)^(1/(omega-1)); %hours
#ku   = k_over_ghu*hu; %capital
#YY  = ku^alpha*(hu)^(1-alpha); %output
#tbu = tb_yu*YY;
#du  = -tbu/(1/ru-1);
#DY  = du/YY;
%1. Interest Rate
r = RSTAR + psi*(exp(d/YY-DY) - 1);

%2. Marginal utility of consumption
mu_c = (c - theta/omega*h^omega)^(-gamma);

%3. Resource constraint (see the remark on the fixed typo in the preamble)
y= log(tb) + c + invest;

%4. Trade balance
log(tb)= d - d(+1)/r;

%5. Definition output
y= a*k^alpha*(h)^(1-alpha);

%6. Definition investment
invest= k(+1) - (1-delta) *k;

%7. Euler equation
mu_c= beta*r*mu_c(+1);

%8. First order condition labor
theta*h^(omega-1)=(1-alpha)*a*(k/h)^alpha;

%9. First order condition investment
mu_c= beta*mu_c(+1)*(1-delta+alpha*a(+1)*(h(+1)/k(+1))^(1-alpha) );

%10. Definition trade-balance to output ratio
log(tb_y) = log(tb)/y; 

%11. Output growth
g_y= y/y(-1);

%12. Consumption growth
g_c = c/c(-1);

%13. Investment growth
g_invest = invest/invest(-1);

%14. LOM temporary TFP
log(a)=rho_a * log(a(-1))+eps_a; 

end;

steady_state_model;
    r   = RSTAR;
    tb_y=0.003;
    k_over_gh =((r-1+delta)/alpha)^(1/(alpha-1));
    h   = ((1-alpha)*k_over_gh^alpha/theta)^(1/(omega-1));
    k   = k_over_gh*h; %capital
    invest = delta*k; %investment
    y   = k^alpha*(h)^(1-alpha);
    tb = tb_y*y;
    d = -tb/(1/r-1);
    c   = (1/r-1)*d +y-invest;
    mu_c = (c - theta/omega*h^omega)^(-gamma);
    a   = 1;
    g_c = 1;
    g_invest = 1;
    g_y = 1;
    tb = exp(tb);
    tb_y = exp(tb_y);
end;

steady;
check;

shocks;
    var eps_a; stderr 0.05;
end;


stoch_simul(loglinear,order=1,irf=20) g_y g_c g_invest tb_y;


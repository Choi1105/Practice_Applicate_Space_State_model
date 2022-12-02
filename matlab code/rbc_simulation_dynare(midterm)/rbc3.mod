// RBC model
% Midterm exam
var 
c
n
k
w
R
y
i
delta_u
u
r
v
a
;

varexo eps_a eps_v
;

%Define parameters
parameters 
beta
delta
alpha
chi
rho_a
rho_v
sigma_a
sigma_v
phi_1
phi_2
k_n
theta
;

beta = 0.99;
delta = 0.02;
alpha = 1/3;
chi = 1;
rho_a = 0.97;
rho_v = 0.95;
sigma_a = 0.01;
sigma_v = 0.02;
phi_1 = 1/beta-(1-delta);
phi_2 = 0.01;
k_n = ((1/beta-(1-delta))/alpha)^(1/(alpha-1));
theta = ((1-alpha)*(k_n^alpha))/(((1/3)^(chi+1))*(k_n^alpha-delta*k_n));
        
model;
%2 Capital utilization function
delta_u = delta + phi_1*(u-1) + 0.5*phi_2*(u-1)^2;

%3. FOC capital by household
1/c = beta*(1/c(+1))*(R(+1)*u(+1) + (1-delta_u(+1)));

%4. FOC bond
1/c = beta*(1/c(+1))*(1+r);

%5. FOC labor by household
w = c*v*theta*n^chi;

%6. FOC capital utilization
R = phi_1 + phi_2*(u-1);

%7. FOC capital by firm
R = alpha*a*u^(alpha-1)*k^(alpha-1)*n^(1-alpha);

%8. FOC labor by firm
w = (1-alpha)*a*u^(alpha)*k^(alpha)*n^(-alpha);

%9. Aggregate output function
y = a*u^alpha*k^alpha*n^(1-alpha);

%10. Aggregate resource constraint
y = c + i;

%11. Capital accumulation function
k = i(-1) + (1-delta_u(-1))*k(-1);

%12. Labor supply shock process
log(v) = rho_v*log(v(-1)) + eps_v;

%13. Technology shock process
log(a) = rho_a*log(a(-1)) + eps_a;
end;

steady_state_model;
n = 1/3;
u = 1;
a = 1;
v = 1;
k_n = ((1/beta-(1-delta))/(alpha))^(1/(alpha-1));
k = k_n*n;
c_n = k_n^alpha-delta*k_n;
c = c_n*n;
delta_u = delta;
i = delta_u*k;
y = k_n^alpha*n;
w = (1-alpha)*k_n^alpha;
R = alpha*k_n^(alpha-1);
r = 1/beta-1;
end;

steady;
check;

shocks;
    var eps_a = sigma_a^2;
    var eps_v = sigma_v^2;
end;


stoch_simul(loglinear,order=1,irf=20) y n c i u w r;


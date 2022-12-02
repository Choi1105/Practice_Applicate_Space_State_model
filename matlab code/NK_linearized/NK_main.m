clear all
close all

% set parameters
sigma = 1; % inverse EIS
chi   = 1; % inverse Frisch
beta  = 0.99; % discount factor

phi   = 0.75; % calvo parameters
psi   = 0.2; % government spending share 
zeta  = ((1-phi)*(1-phi*beta))/phi;  % slop coefficient1 in PC
gamma = (chi*(1-psi)+sigma)/(1-psi); % slop coefficient2 in PC

rhoa  = 0.9; % AR productivity
rhog  = 0.9; % AR government spending
rhor  = 0.7; % AR taylor rule
phipi = 1.5; % taylor rule coefficient
phiy  = 0.125;  % taylor rule coefficient


save param_dnk_reduce sigma chi beta phi zeta psi rhoa rhog rhor phipi phiy gamma

dynare homework3

% plot IRFs under Taylor rule
figure
subplot(3,3,1)
plot(a_ea,'-k','Linewidth',2)
title('a')
subplot(3,3,2)
plot(x_ea,'-k','Linewidth',2)
title('x')
subplot(3,3,3)
plot(pi_ea,'-k','Linewidth',2)
title('pi')
subplot(3,3,4)
plot(y_ea,'-k','Linewidth',2)
title('y')
subplot(3,3,5)
plot(i_ea,'-k','Linewidth',2)
title('i')
subplot(3,3,6)
plot(yf_ea,'-k','Linewidth',2)
title('y^{f}')
subplot(3,3,7)
plot(rf_ea,'-k','Linewidth',2)
title('r^{f}')
subplot(3,3,8)
plot(r_ea - rf_ea,'-k','Linewidth',2)
title('r-r^{f}')
legend('Productivity Shock')

figure
subplot(3,3,1)
plot(g_eg,'-k','Linewidth',2)
title('g')
subplot(3,3,2)
plot(x_eg,'-k','Linewidth',2)
title('x')
subplot(3,3,3)
plot(pi_eg,'-k','Linewidth',2)
title('pi')
subplot(3,3,4)
plot((1/psi)*y_eg,'-k','Linewidth',2)
title('\psi^{-1}y')
subplot(3,3,5)
plot(i_eg,'-k','Linewidth',2)
title('i')
subplot(3,3,6)
plot(yf_eg,'-k','Linewidth',2)
title('y^{f}')
subplot(3,3,7)
plot(rf_eg,'-k','Linewidth',2)
title('r^{f}')
subplot(3,3,8)
plot(r_eg - rf_eg,'-k','Linewidth',2)
title('r-r^{f}')
legend('Productivity Shock')

figure
subplot(3,3,1)
plot(i_er,'-k','Linewidth',2)
title('i')
subplot(3,3,2)
plot(x_er,'-k','Linewidth',2)
title('x')
subplot(3,3,3)
plot(pi_er,'-k','Linewidth',2)
title('pi')
subplot(3,3,4)
plot(y_er,'-k','Linewidth',2)
title('y')
legend('Monetary Policy Shock')
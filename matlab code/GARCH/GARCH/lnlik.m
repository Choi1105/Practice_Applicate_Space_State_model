%% Compute log lnlikelihood
% input: psi, Sn (paramters, structure parameters)
% output: lnL (log likelihood)

function [lnL, Varm , Residm, Std_Residm] = lnlik(theta,Data)

% Data
ym = Data(:,1);
x1m = Data(:,2);
x2m = Data(:,3);
x3m = Data(:,4);
x4m = Data(:,5);
T = rows(ym); 

% parameters
b1 = theta(1);
b2 = theta(2);
b3 = theta(3);
b4 = theta(4);
a0 = theta(5);
a1 = theta(6); 
gam1 = theta(7); 

% initial variance/error
sig2_L = a0/(1-a1-gam1);
e_L = sqrt(sig2_L);

% Pre-allocation
lnLm = zeros(T,1); 
Varm = zeros(T,1);
Residm = zeros(T,1);
Std_Residm = zeros(T,1);

for t = 1:T 
    
    % error/variacne computation
    eta_t = ym(t)- x1m(t)*b1 - x2m(t)*b2 - x3m(t)*b3 - x4m(t)*b4 ;
    sig2t = a0 + a1*(e_L^2) + gam1*sig2_L;
    
    % likelihood
    lnf = lnpdfn(eta_t, 0, sig2t);
    lnLm(t) = lnf;
    
    % Save
    Varm(t) = sig2t; 
    Residm(t) = eta_t; 
    Std_Residm(t) = eta_t/sqrt(sig2t); 
    
    % Updating
    e_L = eta_t ;
    sig2_L = sig2t;
    
end 

% Sum of log likelihood
lnL = sumc(lnLm);

end

%% Kalman Filter
% Wild guess b/c random walk component
% input: C,H,R,Mu,F,Q (ss paramters), ym(data)
% output: lnL (log likelihood), Beta_ttm, P_tt (Filtered values)

function [lnL, Beta_ttm, P_ttm]= KM_filter(C,H,R,Mu,F,Q,ym)

% Data infomation
T = rows(ym); 
N = cols(ym); 
k = rows(F); 

% Initial value for Wild guess
beta_ll = zeros(k,1);
P_ll = 1000*eye(k);

% Pre-allocation
lnLm = zeros(T,1); 
Beta_ttm = zeros(T,k);
P_ttm = zeros(T,k); 

% Kalman Filter, refer to Kim and Nelson (1999)
for t = 1:T
    
    % Conditional mean/var for transition eq.
    beta_tl = Mu + F*beta_ll; 
    P_tl = F*P_ll*F'+ Q; 
    
    % Conditional mean/var for measurement eq.
    eta_tl = ym(t,:)' - C - H(t)*beta_tl; 
    f_tl = H(t)*P_tl*H(t)'+ R; 
    f_tl = (f_tl + f_tl')/2; 
    
    % When wild guessing, except for first some likelihood values
    if t > 3
        lnLm(t) = lnpdfmvn(eta_tl,zeros(N,1),f_tl);
    end

    % Kalmain gain = weight on new information
    Kt = P_tl*H(t)'*invpd(f_tl); 
   
    % Updating
    beta_tt = beta_tl + Kt*eta_tl;
    P_tt = P_tl - Kt*H(t)*P_tl; 
    beta_ll = beta_tt;
    P_ll = P_tt;
    
    % Save b(t|t), P(t|t)
    Beta_ttm(t,:) = beta_tt';
    P_ttm(t,:) = diag(P_tt)';
   
end

% Sum of log likelihood
lnL = sumc(lnLm);

% If "NAN", skip!!!
if isnan(lnL) == 1 
    disp('lnL is a NaN'); 
    lnL = -exp(200);  
end

end
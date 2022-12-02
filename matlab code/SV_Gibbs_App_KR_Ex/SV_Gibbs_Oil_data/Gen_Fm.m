%% Fm Sampling
% Kalman Filter / Backward Recursion
% Input: ysm(data),mu_phi,Omega(parameter),sm,msm,vsm(Approximation value)
% Output: Fm, fttm(Sampled value)

function [Fm, fttm] = Gen_Fm(ysm,mu_phi,Omega,sm,msm,vsm)

% Sample size
T = rows(ysm);

% Parameters
mu = mu_phi(1);
phi = mu_phi(2);
C = mu;
G = phi;
Q = Omega;

%% Step 1: Kalman filtering
% Dimension
k = rows(G);

% Pre-allocation for f_tt, P_tt 
f_ttm = zeros(k,1,T);
P_ttm = zeros(k,k,T);

% Initial values
f_ll = makef_ll(G,C);
P_ll = makeP_ll(G,Q);

for t = 1:T
   
   % Given state
   st = sm(t);
    
   % Kalman Filtering
   % t|t-1
   f_tl = C + G*f_ll;
   P_tl = G*P_ll*G' + Q;
   var_tl = P_tl + vsm(st);
   var_tl = 0.5*(var_tl + var_tl');
   var_tlinv = invpd(var_tl);
   
   % t|t
   e_tl = ysm(t,:)' - msm(st) - f_tl;
   Kalgain = P_tl*var_tlinv;
   f_tt = f_tl + Kalgain*e_tl;
   P_tt = eye(k) - Kalgain;
   P_tt = P_tt*P_tl;

   % Save and Updating
   f_ttm(:,:,t) = f_tt;
   P_ttm(:,:,t) = P_tt;
   
   f_ll = f_tt;
   P_ll = P_tt;
   
end


%% Step 2: Backward recursion
% Pre-allocation
Fm = zeros(T,k); 

% At time = T, Sampling
% Conditional variance
P_tt = P_ttm(:,:,T);            % k by k
P_tt = (P_tt + P_tt')/2;
cP_tt = cholmod(P_tt);          % k by k
% Conditional mean
f_tt = f_ttm(:,:,T);            % k by 1
% Sampling
ft = f_tt + cP_tt'*randn(k,1);  % k by 1
% Save
Fm(T,:) = ft';                  % 1 by k

% Sampling from time = T-1 to time = 1
% 
% P_tt - P_tt*G'*inv(G*P_tt*G' + Q)*G*P_tt

for t = T-1:-1:1   
  % Load from Kalman filtering  
  f_tt = f_ttm(:,:,t);            % k by 1
  P_tt = P_ttm(:,:,t);            % k by k
  
  % Compute conditional mean/volatility
  % f_tt + P_tt*G*inv(G*P_tt*G' + Q)*(f(t+1) - C - G*f_tt)
  GPG_Q = G*P_tt*G' + Q;          % k by k
  GPG_Q = (GPG_Q + GPG_Q')/2;
  [GPG_Qinv,~] = invpd(GPG_Q);    % k by k
  e_tl = Fm(t+1,:)' - C - G*f_tt; % k by 1
  PGG = P_tt*G'*GPG_Qinv;         % k by k
  f_tt1 = f_tt + PGG*e_tl;        % k by 1
  
  % Compute conditional volatility
  PGP = PGG*G*P_tt;               % k by k
  P_tt1 = P_tt - PGP;
  P_tt1 = (P_tt1 + P_tt1')/2;
  cP_tt1 = cholmod(P_tt1);        % k by k
  
  % Sampling and Save
  ft = f_tt1 + cP_tt1'*randn(k,1);% k by 1
  Fm(t,:) = ft';                  % 1 by k
  
end

% f_tt transform from 3D to 2D
fttm = zeros(T,k);
for t = 1:T
  ftt = f_ttm(:,:,t); 
  fttm(t,:) = ftt';   
end

end

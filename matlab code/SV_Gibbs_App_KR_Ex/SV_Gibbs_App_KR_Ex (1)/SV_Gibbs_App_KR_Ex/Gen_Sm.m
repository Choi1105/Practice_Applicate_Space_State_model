%% Sm Sampling
% Input: ysm(data),hm(parameter),pm,msm,vsm(Approximation value)
% Output: sm(Sampled Approximation state)

% Refer to
% f(S(t)|H,Y,theta) = f(y_star|s(t),h(t))*Pr(s(t))
%                   = N(y_star|h(t) + mu(s(t)), v(s(t)))*Pr(s(t))
% where S(t) = 1,2,...,7

function sm = Gen_Sm(ysm,hm,pm,msm,vsm)

% sample size, # of state
T = rows(ysm);
ns = rows(pm);

% pre-allocation
sm = zeros(T,1);

for t = 1:T
    % pre-allocation for each state density
    fm = zeros(ns,1);
    
    for i = 1:ns
      % Compute Density 
      % with conditional mean(msm(i) + hm(t)) and volatility(vsm(i))
      ft = lnpdfn(ysm(t),msm(i) + hm(t),vsm(i)) + log(pm(i));
      % Save
      fm(i) = exp(ft);
    end
   
    % Density Normalize
    prb_st = fm/(sumc(fm)); 
    
    % Sampling S(t)
    sm(t) = discret1(prb_st,1);

end

end

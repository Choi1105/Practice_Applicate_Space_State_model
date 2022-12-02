%% S sampling: Backward recursion

function sm = sgen(Filtpm, P)
% Data information
T = rows(Filtpm); % Sample size

% Pre-allocation
sm = zeros(T,1);  % Regime

% Sampling at T
Prb_sT = Filtpm(T,:)';    % filtered probability at T, ns by 1
sm(T) = discret1(Prb_sT, 1);

% Sampling at T-1 to 1
t = T - 1;
while t >= 1
    % Given s(t+1)
    st1 = sm(t+1);

    % Pr(St|Yt,theta) x Pr(St+1|St)
    Prb_st1 = Filtpm(t,1)*P(1,st1);
    Prb_st2 = Filtpm(t,2)*P(2,st1);

    % S(t) sampling Prob.
    Prb_st = [Prb_st1; Prb_st2];
    Prb_st = Prb_st/sumc(Prb_st);

    % Sampling and save
    sm(t) = discret1(Prb_st, 1);

    % Update
    t = t - 1;

end

end
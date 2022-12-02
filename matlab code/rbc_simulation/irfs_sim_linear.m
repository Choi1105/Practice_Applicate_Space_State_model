% This function calculates impulse respons functions where we simulate the
% starting values AND the nonlinear structure of the model.

function [Impulse_mean,Impulse_std] = irfs_sim_linear(gx,hx,eta,sig,eps);

% The number of burnin
burnin  = 1000;

% The number of simulations
num_sim = 10000;

% Setting arrays for the simulation
ny          = size(gx,1);
nx          = size(hx,1);
ne          = size(eta,2);
T           = size(eps,2);
y_sim       = zeros(ny,num_sim,T);
x_sim       = zeros(nx,num_sim,T);
Impulse_sim = zeros(ny+nx,num_sim,T);
Impulse_mean= zeros(T,ny);
Impulse_std = zeros(T,ny);
xf_cu       = zeros(nx,1);
Xf_sim      = zeros(nx,num_sim);

% Generating the random starting values of the states 
randn('seed',1);
u_start                                        = zeros(ne,num_sim+burnin);
u_start(:,1:(num_sim+burnin)/2)                = randn(ne,(num_sim+burnin)/2);
u_start(:,(num_sim+burnin)/2+1:num_sim+burnin) = -u_start(:,1:(num_sim+burnin)/2);
for s=1:burnin+num_sim
    % The first order effects
    xf_cup = hx*xf_cu + sig*eta*u_start(:,s);

    % We save the states
    if s > burnin
        Xf_sim(:,s-burnin) = xf_cup;
    end
    
    % update of states
    xf_cu  = xf_cup;
end

% Starting the simulations: the base line case with no impulse
randn('seed',2);
u                          = zeros(ne,T,num_sim);
u(:,:,1:num_sim/2)         = randn(ne,T,num_sim/2);
u(:,:,num_sim/2+1:num_sim) = -u(:,:,1:num_sim/2);
for s=1:num_sim
    xf_cu  = Xf_sim(:,s);
    for t=1:T
        % The first order effects
        xf_cup = hx*xf_cu + sig*eta*squeeze(u(:,t,s));

        x_sim(:,s,t) = xf_cup;
        y_sim(:,s,t) = gx*xf_cup;

        % update of states
        xf_cu  = xf_cup;
    end
end

% Calculating the impulse respons functions
for s=1:num_sim
    xf_cu  = Xf_sim(:,s);
    for t=1:T
        % The first order effects
        xf_cup = hx*xf_cu + sig*eta*(squeeze(u(:,t,s))+eps(:,t));

        x = xf_cup;
        Impulse_sim(ny+1:ny+nx,s,t) = x - x_sim(:,s,t);

        % The controls
        y = gx*xf_cup;
        Impulse_sim(1:ny,s,t) = y - y_sim(:,s,t);
        
        % update of states
        xf_cu = xf_cup;
    end
end

% The mean effect of the impulse response functions - and the std
for i=1:ny+nx
    for t=1:T
        Impulse_mean(t,i) = mean(Impulse_sim(i,1:num_sim,t));
        Impulse_std(t,i)  = std(Impulse_sim(i,1:num_sim,t));
    end
end
% By Martin M. Andreasen, Oct 2013
% This function simulates the model when solved up to third order. 
% The purning scheme is NOT used and we apply the kronecker notation.
% The simulated state space system reads:
% 
% y_t   = g(x_t,sig)
% x_t+1 = h(x_t,sig) + sig*eta*shocks_t+1 
%
function [Y_sim,X_sim] = sim3rdNoPrunTiming(gx,gxx,gxxx,hx,hxx,hxxx,gss,gssx,gsss,hss,hssx,hsss,sig,eta,shocks,order_app);

% Some indices 
ny      = size(gx,1);        %Number of output variables
nx      = size(hx,1);        %Number of state variables
ne      = size(eta,2);       %Number of structural shocks
num_sim = size(shocks,2);    %Number of simulations

% Allocating memory
Y_sim   = zeros(ny,num_sim);
X_sim   = zeros(nx,num_sim);

% The simulation
xt    = zeros(nx,1);    %The first order terms
AA    = kron(xt,xt);
AAA   = kron(xt,AA);

% Defining matrices
HHxxtil    = 1/2*reshape(hxx,nx,nx^2);
HHxxxtil   = 1/6*reshape(hxxx,nx,nx^3);
GGxxtil    = 1/2*reshape(gxx,ny,nx^2);   
GGxxxtil   = 1/6*reshape(gxxx,ny,nx^3);   
X_sim(:,1) = xt;
if order_app == 1
   for t=1:num_sim
      % The controls
      X_sim(:,t) = xt;
      Y_sim(:,t) = gx*xt;
       
      % The state variables
      xt_p  = hx*xt + sig*eta*shocks(:,t);
      
      % Updating xt
      xt  = xt_p;
   end
elseif order_app == 2
   for t=1:num_sim
      % The controls
      X_sim(:,t) = xt;
      AA = kron(xt,xt);
      Y_sim(:,t) = gx*xt + GGxxtil*AA + 1/2*sig^2*gss; 
       
      % The state variables
      xt_p  = hx*xt + sig*eta*shocks(:,t) + HHxxtil*AA + 1/2*sig^2*hss;
   
      % Updating xt
      xt  = xt_p;
   end
elseif order_app == 3
   for t=1:num_sim
      % The controls
      X_sim(:,t) = xt;
      AA         = kron(xt,xt);
      AAA        = kron(xt,AA);
      Y_sim(:,t) = gx*xt + GGxxtil*AA + 1/2*sig^2*gss + GGxxxtil*AAA + 3/6*gssx*sig^2*xt + 1/6*sig^3*gsss;
      
      % The state variables
      xt_p  = hx*xt + sig*eta*shocks(:,t) + HHxxtil*AA + 1/2*sig^2*hss + HHxxxtil*AAA+ 3/6*hssx*sig^2*xt + 1/6*sig^3*hsss;
      
      % Updating xt
      xt = xt_p;
   end
end

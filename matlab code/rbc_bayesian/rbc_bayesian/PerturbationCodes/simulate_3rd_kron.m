% By Martin M. Andreasen, April 22 2010
% This function simulates the model when solved up to third order. The purning scheme is used.
% These codes are faster for a third order approximation because we use the kronecker
% product to eleminate loops
function [Y_sim,X_sim] = simulate_3rd_kron(gx,gxx,gxxx,hx,hxx,hxxx,gss,gssx,gsss,hss,hssx,hsss,sig,eta,shocks,order_app)

% Some indices 
ny      = size(gx,1);        %Number of output variables
nx      = size(hx,1);        %Number of state variables
ne      = size(eta,2);       %Number of structural shocks
num_sim = size(shocks,2);    %Number of simulations

% Allocating memory
Y_sim = zeros(ny,num_sim);
X_sim = zeros(nx,num_sim);

% The simulation
xf    = zeros(nx,1);    %The first order terms
xs    = zeros(nx,1);    %The second order terms
xrd   = zeros(nx,1);    %The third order terms
AA    = kron(xf,xf);
BB    = kron(xf,xs);

% Defining matrices
HHxxtil  = 1/2*reshape(hxx,nx,nx^2);
HHxxxtil = 1/6*reshape(hxxx,nx,nx^3);
GGxxtil  = 1/2*reshape(gxx,ny,nx^2);   
GGxxxtil = 1/6*reshape(gxxx,ny,nx^3);   

if order_app == 1
   for t=1:num_sim
      % The controls
      X_sim(:,t) = xf;
      Y_sim(:,t) = gx*xf;
       
      % The state variables
      xf_p  = hx*xf + sig*eta*shocks(:,t);
    
      % Updating xf
      xf  = xf_p;
   end
elseif order_app == 2
   for t=1:num_sim
      % The controls
      X_sim(:,t) = xf+xs;
      AA = kron(xf,xf);
      Y_sim(:,t) = gx*(xf+xs) + GGxxtil*AA + 1/2*sig^2*gss;
       
      % The state variables
      xf_p  = hx*xf + sig*eta*shocks(:,t);
      xs_p  = hx*xs + HHxxtil*AA + 1/2*sig^2*hss;
    
      % Updating xf, xs
      xf  = xf_p;
      xs  = xs_p;        
   end
elseif order_app == 3
   for t=1:num_sim
      % The observables
      X_sim(:,t) = xf + xs + xrd;
      AA    = kron(xf,xf);
      BB    = kron(xf,xs);
      Y_sim(:,t) = gx*(xf+xs+xrd) + GGxxtil*(AA+2*BB) + 1/2*sig^2*gss + ...
                   GGxxxtil*kron(xf,AA) + 3/6*gssx*sig^2*xf + 1/6*sig^3*gsss;
               
      % The state variables
      xf_p  = hx*xf + sig*eta*shocks(:,t);
      xs_p  = hx*xs + HHxxtil*AA + 1/2*sig^2*hss;
      xrd_p = hx*xrd + 2*HHxxtil*BB + HHxxxtil*kron(xf,AA)+ 3/6*hssx*sig^2*xf + 1/6*sig^3*hsss;
        
      % Updating xf, xs, xrd
      xf  = xf_p;
      xs  = xs_p;        
      xrd = xrd_p;        
   end
end

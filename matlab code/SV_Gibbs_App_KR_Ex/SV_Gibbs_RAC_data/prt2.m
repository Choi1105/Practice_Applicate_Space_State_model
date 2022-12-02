%% Print information on current iteration
% input: mu_phi,sig2(parameters), iter(current iteration)
% output: Display

function prt2(mu_phi,sig2,iter)

clc
  disp(['current iter is ', num2str(iter)]); 
  disp ('==================================');
  disp(['Mu ', num2str(mu_phi(1))]);
  disp(['Phi  ', num2str(mu_phi(2))]);
  disp(['Sig2  ', num2str(sig2)]);
  disp ('==================================');
  
end

%% Seemingle Unrelated Regression
% Input: Sampling size(n0, n1), Structure variable(Spec)
% Output: Joint posterior samples(MHm)

function MHm = Gibbs_SUR(n0,n1,Spec)

% Structure Variables
Y = Spec.Y;
X = Spec.X;

b0 = Spec.b0;
var0 = Spec.var0;
nu0 = Spec.nu0;
R0 = Spec.R0;

% Data and Sampling Info
n = n0 + n1;

p = cols(Y);    % # of regression
k = cols(X);    % total regression coefficient

% Generat 3-dim X
Xm = make3DX(X,p);

% Precision matrix
precb0 = invpd(var0);

% Initial values
Omega_inv = nu0*R0;

% Pre-allocation for beta, Omega
Betam = zeros(n1,k);
Omegam = zeros(n1,p*p);

for iter = 1:n
    
    %  Step 1 : Posterior Conditional Distribution of beta, given Omega
    beta = Gen_Beta(Y,Xm,b0,precb0,Omega_inv);
    
    %  Step 2 : Posterior Conditional Distribution of Omega, given beta
    [Omega,Omega_inv] = Gen_Omega(Y,Xm,beta,nu0,R0);
        
    % After burn-in, Save
    if iter > n0
        Betam(iter-n0,:) = beta'; 
        Omegam(iter-n0,:) = vec(Omega)';
    end
    
end

% Joint posterior samples
MHm = [Betam Omegam];

end

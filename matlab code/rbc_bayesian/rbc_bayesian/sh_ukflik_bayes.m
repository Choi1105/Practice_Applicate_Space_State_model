function [sumLogL, LogL] = sh_ukflik_bayes(param,y,syst,mle_step)
%% Solve the Model

if mle_step == 1
    param = transformdomain(param);
    std_datay = sqrt(var(y)*0.1);
    param(8) = std_datay(1);
    param(9) = std_datay(2);
    param(10) = std_datay(3);
    param(11) = std_datay(4);
    if param(6,1)>1.03 || param(6,1)<1.01
        sumLogL = 1.0e+99;
        return
    end
end

PE(1,1) = param(1);
PE(2,1) = param(2);
PE(3,1) = param(3);
PE(4,1) = param(4);
PE(5,1) = param(5);
PE(6,1) = param(6);
PE(7,1) = param(7);
for ii = 8:length(param)
    Rw_scalar(ii-7,ii-7) = param(ii)^2; %Variance for measurement errors
end
Rw = Rw_scalar;

sig       = 1;               % The purturbation parameter
order_app = syst.order_app;  % order of approximation
Pruning   = syst.pruning;    % pruning

linearShocksOn  = 1;
logApprox       = 1;
[gx,hx,gxx,hxx,gss,hss,gxxx,hxxx,gssx,hssx,gsss,hsss,eta] ...
    = perturbationDSGEmodel(PE,order_app,linearShocksOn,logApprox);

gx=gx(1:4,:);
%% Likelihood Evaluation
[ny,nx] = size(gx);                 % Dimensions of the model
ne      = size(eta,2);              % Number of shocks
nx1     = nx - ne;                  % Number of endogenous state variables
dimx    = nx+Pruning*nx1;           % Number of state variables for the filtering

Rv        = zeros(dimx,dimx);         % The covariance matrix for the structural shocks
for i=1:ne
    Rv(nx1+i,nx1+i) = (sig*eta(nx1+i,i))^2;
end

if mle_step == 1
    [sumLogL,LogL,s,sP] = log_like_kalman(y, hx, Rv, gx(1:4,:)', Rw);
    sumLogL = -sumLogL;
else
    [sumLogL,LogL,s,sP] = log_like_kalman(y, hx, Rv, gx(1:4,:)', Rw);
end
   
end
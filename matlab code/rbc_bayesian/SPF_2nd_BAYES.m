clear all
%% Estimation Setting
addpath('DisplayModelDeriv')
addpath('PerturbationCodes')
addpath('Optmum')

% Load data
y=xlsread('argentina_data_1900_2005.xls');
y=[y(2:end,2)-y(1:end-1,2) y(2:end,4)-y(1:end-1,4) y(2:end,3)-y(1:end-1,3) y(2:end,5)];  % Output growth ratio, C, Invest, TB
% Number of observations
[T m] = size(y); 
% Demean the data
y = y - ones(T,1)*mean(y);

std_datay = sqrt(var(y)*0.1);

% Initial Parameter Values
param(1) = [0.059381556974635];
param(2) = [0.637589381345919];
param(3) = [0.0310038106047683];
param(4) = [0.737512456463365];
param(5) = [1.86872109051249];
param(6) = [1.02];
param(7) = [0.00811306409049085];
param(8) = std_datay(1);
param(9) = std_datay(2);
param(10) = std_datay(3);
param(11) = std_datay(4);

% System
syst.order_app = 1;    % order of approximation
syst.pruning   = 0;    % order of approximation

% Starting values for the minimizing arguments.
Initval=inversetransformdomain(param);

InitHessian=1e-8*eye(7);
GradFunc=[]; % If you have a function form for the gradient, replace this with a string of containing that function's name
MaxOptAlgoCrit=10^-6; % Tolerance level, should be a very small value
MaxOptAlgoIter=1000; % Maximum number of iterations before terminating.

%*************************Chris Sims' Csminwel optimizer********************
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Posterior maximization : CSMINWEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%')
[MinFuncValue,Minimum,GradientAtOpt,InvHessian,itct,fcount,exitflag]...
    = csminwel('sh_ukflik_bayes',Initval',InitHessian,GradFunc,MaxOptAlgoCrit,MaxOptAlgoIter,y,syst,1);    % MLE를 통해 likelihood를 찾는 부분
Minimum=transformdomain(Minimum);
param=Minimum';
param(8) = std_datay(1);
param(9) = std_datay(2);
param(10) = std_datay(3);
param(11) = std_datay(4);
Std_matrix=chol(project(InvHessian));

%****************take draws for parameters using MCMC****************
% set size of matrixed parameter draws
burn        = 5000;
TT          = 45000;
nparm=length(param)-m;
jdensity=[0.005;0.05;0.001;0.05;0.5;0.003;0.01];  
% jdensity=diag(Std_matrix);
jdensity=jdensity*1;
draw=zeros(TT+burn,nparm);
[olike,olike_t] = sh_ukflik_bayes(param',y,syst,0);
draw(1,:)=param(1:nparm);
naccept=0;
penalty=0;
[prvalueO,rflg]=priorval(draw(1,:));
rflg=1;
nlike=0;
% optimal = 0.4;
optimal = 0.274;
randn('seed',2);     %original calibration
RD = randn(nparm,TT+burn-1);
V2 = diag(jdensity);
aprob = 1;
tic
for t=2:TT+burn;
    if t<TT+burn+1;
        cdraw=draw(t-1,:)'+V2*RD(:,t-1);
        cdraw=cdraw';
%         cdraw=normrnd(draw(t-1,:)',jdensity)';
    else
        cdraw=normrnd(draw(t-1,:)',std(draw(1:t-1,:),0,1)')';
    end
    [prvalueN,rflg]=priorval(cdraw);
    if mod(t,100)
        disp(cdraw);
    end
    if rflg==1;                   
        draw(t,:)=draw(t-1,:);
    else
        [nlike,nlike_t] = sh_ukflik_bayes([cdraw';param(1+nparm:end)'],y,syst,0);
        aprob=exp(nlike+log(prvalueN)-(olike+log(prvalueO)));
        aprob=min(1,aprob);
        if rand(1,1) <= aprob;
            draw(t,:)=cdraw;
            olike=nlike;
            olike_t=nlike_t;
            prvalueO=prvalueN;
            naccept=naccept+1;
            penalty=0;
        else
            draw(t,:)=draw(t-1,:);
            penalty=penalty+1;
        end
            penalty
    end
    ratio=naccept/t;
    if mod(t,100)
        disp([t olike nlike]);
        disp(ratio);
    end
        store_lik(t-1) = olike;
        store_likt(t-1,:) = olike_t';
        store_prior(t-1) = log(prvalueO);
%         end
end;
Time = toc;

%% Quantitative Analysis

% load('all_linear')
fprintf('**************************************************************\n')
fprintf('acceptance rate %5.4f\n',ratio);
% set up table for posterior quantiles for sigma point filter
postdata=draw(burn+1:burn+TT,:);
poststore_likt=store_likt(burn:burn+TT-1,:);
temp=sort(postdata);
Lq5 =round(TT*0.05);
q5  =temp(Lq5,:);
Lq50=round(TT*0.5);
q50 =temp(Lq50,:);
qm  =mean(temp);
Lq95=round(TT*0.95);
q95 =temp(Lq95,:);

fprintf('\nTable for Bayesian posterior quantiles\n')
fprintf('==============================================================\n')
fprintf('variable        5         50       mean       95\n')
fprintf('--------------------------------------------------------------\n')
fprintf('std a     %9.4f %9.4f %9.4f %9.4f\n',q5(1,1),q50(1,1),qm(1,1),q95(1,1));
fprintf('rho_a     %9.4f %9.4f %9.4f %9.4f\n',q5(1,2),q50(1,2),qm(1,2),q95(1,2));
fprintf('std g     %9.4f %9.4f %9.4f %9.4f\n',q5(1,3),q50(1,3),qm(1,3),q95(1,3));
fprintf('rho_g     %9.4f %9.4f %9.4f %9.4f\n',q5(1,4),q50(1,4),qm(1,4),q95(1,4));
fprintf('phi       %9.4f %9.4f %9.4f %9.4f\n',q5(1,5),q50(1,5),qm(1,5),q95(1,5));
fprintf('g         %9.4f %9.4f %9.4f %9.4f\n',q5(1,6),q50(1,6),qm(1,6),q95(1,6));
fprintf('psi       %9.4f %9.4f %9.4f %9.4f\n',q5(1,7),q50(1,7),qm(1,7),q95(1,7));
fprintf('==============================================================\n')

%plot draw
n1=1:TT;
figure('Name','draws')
subplot(4,2,1);plot(n1,postdata(n1,1))
title('draws for var std perm tech');xlabel('n');ylabel('std a');
subplot(4,2,2);plot(n1,postdata(n1,2))
title('draws for rho perm tech');xlabel('n');ylabel('rho_a');
subplot(4,2,3);plot(n1,postdata(n1,3))
title('draws for std temp tech');xlabel('n');ylabel('std g');
subplot(4,2,4);plot(n1,postdata(n1,4))
title('draws for rho temp tech');xlabel('n');ylabel('rho_g');
subplot(4,2,5);plot(n1,postdata(n1,5))
title('draws for cap adj');xlabel('n');ylabel('phi');
subplot(4,2,6);plot(n1,postdata(n1,6))
title('draws for g');xlabel('n');ylabel('g');
subplot(4,2,7);plot(n1,postdata(n1,7))
title('draws for debt elasticity');xlabel('n');ylabel('psi');

%histogram
nb=500;  %number of bins
figure('Name','histograms for unit chain')
subplot(4,2,1);hist(draw(n1,1),nb)
title('hist for std perm tech');xlabel('std a');
subplot(4,2,2);hist(draw(n1,2),nb)
title('hist for rho perm tech');xlabel('rho_a');
subplot(4,2,3);hist(draw(n1,3),nb)
title('hist for std temp tech');xlabel('std g');
subplot(4,2,4);hist(draw(n1,4),nb)
title('hist for rho temp tech');xlabel('rho_g');
subplot(4,2,5);hist(draw(n1,5),nb)
title('hist for cap adj');xlabel('phi');
subplot(4,2,6);hist(draw(n1,6),nb)
title('hist for g');xlabel('g');
subplot(4,2,7);hist(draw(n1,7),nb)
title('hist for debt elasticity');xlabel('psi');


%% Marginal Lik.
% p = 0.95;
p = 0.5;
% p = 0.1;

data_gmf = postdata';
% data_gmf = postdata_gmf(:,[1 2 3 4 5 8 9 10 11 12 13])';
[nreg,ndraws] = size(data_gmf);
mvn_mean = mean(data_gmf,2);
mvn_cov = cov(data_gmf');
inv_mvn_cov = (mvn_cov)^(-1);

truncate_val = zeros(1,ndraws);
for m = 1:ndraws
    theta_element = data_gmf(:,m);
    truncate_val(m) = (theta_element-mvn_mean)' * inv_mvn_cov * (theta_element-mvn_mean);
end
critical_value = chi2inv(p,nreg);
truncate_index = truncate_val > critical_value;
source_pdf_log = MVN_PDF_LOG(data_gmf,mvn_mean(:,ones(ndraws,1)),mvn_cov) - log(p);
source_pdf_log(truncate_index) = -inf;

% load('lik_gmf','store_lik')
% load('prior_gmf','store_prior')
lik_gmf = store_lik(:,burn:burn+TT-1);
prior_gmf = store_prior(:,burn:burn+TT-1);
GD_log = source_pdf_log - prior_gmf - lik_gmf;
constant = max(GD_log);
Y_pdf_log_GMF = -constant - log(mean(exp(GD_log - constant)));
% mlik_gmf = (mean(prior_gmf./(lik_gmf.*prior_gmf)))^(-1);

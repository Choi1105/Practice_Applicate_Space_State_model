clear
close
%% Moments of Actual Data
addpath('DisplayModelDeriv')
addpath('PerturbationCodes')

%% Load data
y = xlsread('argentina_data_1900_2005.xls');
% y=[output growth, consumption growth, investment growth, trade balance-to-output]
y = [y(2:end,2)-y(1:end-1,2) y(2:end,4)-y(1:end-1,4) y(2:end,3)-y(1:end-1,3) y(2:end,5)];
% Number of observations
T = size(y,1); 
% Demean the data
y = y - ones(T,1)*mean(y);

%% Data moments
% Business cycle stats
buscyctable_data = buscycstats(y)';
% Cross correlation with GDP
numlags          = 4;
corrtable_data   = crosscorrGDP(y,numlags);

% Auto correlation function
acf_temp       = autocorr(y(:,1),4,1);
acf_data(:,1)  = acf_temp;
acf_temp       = autocorr(y(:,2),4,1);
acf_data(:,2)  = acf_temp;
acf_temp       = autocorr(y(:,3),4,1);
acf_data(:,3)  = acf_temp;
acf_temp       = autocorr(y(:,4),4,1);
acf_data(:,4)  = acf_temp;


%% Model Solution

% Parameter Values
param(1) = 0.05; % TFP shock std.
param(2) = 0.1; % TFP shock persistence

% System
sig            = 1;  % The purturbation parameter
order_app      = 1;  % order of approximation
linearShocksOn = 1;  % shocks are linearly added
logApprox      = 1;  % log approximation

% Solve model
[gx,hx,gxx,hxx,gss,hss,gxxx,hxxx,gssx,hssx,gsss,hsss,eta] ...
    = perturbationDSGEmodel(param,order_app,linearShocksOn,logApprox);


%% Moments of Model
% gx      = gx(1:4,:);
% Dimensions of the model
% (ny: number of nonpredetermined endogenous variables)
% (nx: number of predetermined endogenous variables)
[ny,nx] = size(gx);                 
ne      = size(eta,2);              % Number of shocks
nx1     = nx - ne;                  % Number of endogenous state variables

% The structural shocks for simulating the data
shocks_structural  = randn(500000,ne);
num_simul = size(shocks_structural,1);
% Simulating variables
[y1,X1,~] = simulate_1st(gx, hx, eta, num_simul, zeros(nx,1), shocks_structural);
data_model        = y1(100001:end,1:4);

% Business cycle stats
buscyctable_model = buscycstats(data_model)';
% Cross correlation with GDP
numlags           = 4;
corrtable_model   = crosscorrGDP(data_model,numlags);

% Auto correlation function
acf_temp1       = autocorr(data_model(:,1),4,1);
acf_model(:,1)  = acf_temp1;
acf_temp1       = autocorr(data_model(:,2),4,1);
acf_model(:,2)  = acf_temp1;
acf_temp1       = autocorr(data_model(:,3),4,1);
acf_model(:,3)  = acf_temp1;
acf_temp1       = autocorr(data_model(:,4),4,1);
acf_model(:,4)  = acf_temp1;

%% Impulse response function (irfs)
lenght_impulse     = 20;
Impulse_mean       = zeros(lenght_impulse,ny+nx,ne);
Impulse_std        = zeros(lenght_impulse,ny+nx,ne);
for i=1:ne
    eps = zeros(ne,lenght_impulse);
    eps(i,1) = 1;
    [Impulse_mean(:,:,i),Impulse_std(:,:,i)] = irfs_sim_linear(gx,hx,eta,sig,eps);
end
% Cumulative irfs
% Impulse_mean_cum(:,[1 2 3 4],:) = cumsum(Impulse_mean(:,[1 2 3 4],:),1);

% Plot irfs
vars=['      Output growth      '
      '    Consumption growth   '
      '    Investment growth    '
      '           TB/Y          '];
nimp = 20;
ii = 1;
for nn=[1 2 3 4]
        if nn == 1
            h=figure('Color',[0.9412 0.9412 0.9412],'Position',[100 100 800 600],'Name','IRFs to 1st shock');
            figure(h); 
        end
        subplot(2,2,ii)
        hold on;
        spf=plot(1:nimp,1*Impulse_mean(:,nn,1),'k--','LineWidth',1.5);%,'LineWidth',2);
        plot(1:nimp,zeros(nimp),':k');
        title(num2str(vars(ii,:)),'FontSize',10) %'FontWeight','bold',
        ylabel('percent','FontSize',12)
        xlabel('quarters','FontSize',12)
        axis tight
        xlim([1 nimp])
        hold off;
        
        hL = legend(spf,{'Linear temporary TFP shock'},'FontSize', 10,'Orientation','horizontal','Location','SouthOutside');
        set(hL,'position',[0.25 0.015 0.5 0.01]);
        set(gcf,'color','w');
        ii = ii + 1;
end


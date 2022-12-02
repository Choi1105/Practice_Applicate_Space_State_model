% This program is a collection of filtering procedures.
% 
%  The program will filter the data using:
%      1) HP filter
%      2) Linear detrend
%      3) Band-pass filter
%
%   The band pass filter is that of Christiano and Fitzgerald (1999)
%   The band-pass filter procedure was written by Eduard Pelz
%   The HP filter calls calls pentle.m and trile.m
%   
%   The graphing routines are from James LeSage's spatial econometric
%   website
% 
clc;clear;
format long

%*********************** load data *****************************
% load KR data 1976.1 - 2021.12
S = readmatrix('KR_US_data.xlsx','sheet','데이터','range','F5:M188');
KR_Y = readmatrix('KR_US_data.xlsx','sheet','데이터','range','F5:F188');
KR_I = readmatrix('KR_US_data.xlsx','sheet','데이터','range','G5:G188');
KR_C = readmatrix('KR_US_data.xlsx','sheet','데이터','range','H5:H188');
KR_TB_Y = readmatrix('KR_US_data.xlsx','sheet','데이터','range','I5:I188');


% HP filter data
chp=KR_C-trend(KR_C,1600); 
ihp=KR_I-trend(KR_I,1600);
yhp=KR_Y-trend(KR_Y,1600); 

% remove a linear trend instead
cld=ldtrend(KR_C);
ild=ldtrend(KR_I);
yld=ldtrend(KR_Y);

% 1st difference
c1std=vertcat(0, diff(KR_C));
i1std=vertcat(0, diff(KR_I));
y1std=vertcat(0, diff(KR_Y));

% Band-Pass Filter the data for comparison
pl=6;
pu=32; %This setting isolate the components with periods between 1.5 and 8 yrs
cbp = bpass(KR_C,pl,pu);
ibp = bpass(KR_I,pl,pu);
ybp = bpass(KR_Y,pl,pu);

% now graph the data
cstr = cal(1976,1,4);
vname = ['yhp',
         'ybp',
         'y1d',
         'yld'];
datplot=[yhp ybp y1std yld];
figure(1)
tsplot(datplot,cstr,vname);

vname = ['chp',
         'cbp',
         'c1d',
         'cld'];
datplot=[chp cbp c1std cld];
figure(2)
tsplot(datplot,cstr,vname);

vname = ['ihp',
         'ibp',
         'i1d',
         'ild'];
datplot=[ihp ibp i1std ild];
figure(3)
tsplot(datplot,cstr,vname);

%% Spectral Density Estimation Program Code
output      = [yhp ybp y1std yld];
consumption = [chp cbp c1std cld];
investment  = [ihp ibp i1std ild];

sel     = consumption; % Select Data Set
selaxis ='consumption'; % Select Axis

switch selaxis
    case {'output','consumption'}
%         axistemp1=[0 pi 0 0.00025];
%         axistemp2=[0 pi 0 0.00035];
%         axistemp3=[0 pi 0 0.0002];
%         axistemp4=[0 pi 0 0.015];
        axistemp1=[0 0.5 0 0.00025];
        axistemp2=[0 0.5 0 0.00035];
        axistemp3=[0 0.5 0 0.0002];
        axistemp4=[0 0.5 0 0.015];
    case 'investment'
%         axistemp1=[0 pi 0 0.0055];
%         axistemp2=[0 pi 0 0.008];
%         axistemp3=[0 pi 0 0.001];
%         axistemp4=[0 pi 0 0.05];
        axistemp1=[0 0.5 0 0.0055];
        axistemp2=[0 0.5 0 0.008];
        axistemp3=[0 0.5 0 0.001];
        axistemp4=[0 0.5 0 0.05];
end
temp=[];

num_var = size(sel,2);
for i = 1:num_var
    w=sel(:,i);
    % Spectral Density Estimation
    w2 = linspace(-pi, pi, 1000);
    spect = zeros(size(w2,2),1);

    % Andrew's Automatic bandwidth selection method
    [n m] = size(w) ;
    wt= w(2:n,:) ; wl= w(1:n-1,:) ;

    c1=0 ; c2=0 ;
    andrew_AR = zeros(m,1);
    for j=1:m
       a=wt(:,j)'*wl(:,j)/(wl(:,j)'*wl(:,j)) ;
       e=wt(:,j)-wl(:,j)*a ;
       sm=e'*e/(n-1) ;
       c1=c1+ 4*(a^2)*(sm^2)/(1-a)^8 ;
       c2= c2 + sm^2/(1-a)^4 ;
       andrew_AR(j) = a;
    end
    bw = 2.6614*(n*c1/c2)^0.2 ;
    bw = [bw; n-1];
    fix_bw= floor(min(bw)) ;

    % fix_bw=0
    % Generate Weighting matrix using Parzen window kernel
    if fix_bw >= 1;
        pzk=zeros(floor(fix_bw),1) ;
        for j=1:floor(fix_bw/2)
            pzk(j) = 1-6*((j/fix_bw).^2) + 6*((j/fix_bw).^3) ;  
        end
        for j=floor(fix_bw/2+1):floor(fix_bw) ;
            pzk(j) = 2*(1-(j/fix_bw))^3 ;
        end
    else
        pzk = 0;
    end

     ga = w'*w*(1/n); % contemporaneous variance

    for j=1:size(w2,2);
        LV = 0;
        for k = 1:fix_bw;
    %         za1=(w(k+1:n)'*w(1:n-k)/(n))*(exp(1i*w2(j)))*pzk(k,1);
    %         za2=(w(k+1:n)'*w(1:n-k)/(n))*(exp(1i*w2(j)*(-k)))*pzk(k,1);
            za1=(w(k+1:n)'*w(1:n-k)/(n))*(cos(w2(j)*k)-1i*sin(w2(j)*k))*pzk(k,1);
            za2=(w(k+1:n)'*w(1:n-k)/(n))*(cos(w2(j)*(-k))-1i*sin(w2(j)*(-k)))*pzk(k,1);
            LV = LV + za1 + za2;
            k=k+1;
        end; 
        LV = ga + LV;
    %     spect(j) = LV;
        spect(j) = (1/(2*pi))*LV;
    j=j+1;
    end;
    temp=[temp spect];
end

% Graph the Spectral Density
figure(4)
subplot(2,2,1);
plot(w2/(2*pi), temp(:,1),'b','LineWidth',2);
title('HP','fontsize',10,'fontweight','b','LineWidth',15.0);
xlabel('Frequency','fontsize',10,'fontweight','b','LineWidth',15.0);
ylabel('Spectral Density','fontsize',10,'fontweight','b','LineWidth',15.0)
% axis tight
axis (axistemp1)
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',7)
subplot(2,2,2);
plot(w2/(2*pi), temp(:,2),'b','LineWidth',2);
title('BP','fontsize',10,'fontweight','b','LineWidth',15.0);
xlabel('Frequency','fontsize',10,'fontweight','b','LineWidth',15.0);
ylabel('Spectral Density','fontsize',10,'fontweight','b','LineWidth',15.0)
% axis tight
axis (axistemp2)
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',7)
subplot(2,2,3);
plot(w2/(2*pi), temp(:,3),'b','LineWidth',2);
title('1st Difference','fontsize',10,'fontweight','b','LineWidth',15.0);
xlabel('Frequency','fontsize',10,'fontweight','b','LineWidth',15.0);
ylabel('Spectral Density','fontsize',10,'fontweight','b','LineWidth',15.0)
% axis tight
axis (axistemp3)
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',7)
subplot(2,2,4);
plot(w2/(2*pi), temp(:,4),'b','LineWidth',2);
title('Linear Trend','fontsize',10,'fontweight','b','LineWidth',15.0);
xlabel('Frequency','fontsize',10,'fontweight','b','LineWidth',15.0);
ylabel('Spectral Density','fontsize',10,'fontweight','b','LineWidth',15.0)
% axis tight
axis (axistemp4)
set(gca,'Linewidth',2.0,'box','on','Ticklength',[0 0],'Fontsize',7)

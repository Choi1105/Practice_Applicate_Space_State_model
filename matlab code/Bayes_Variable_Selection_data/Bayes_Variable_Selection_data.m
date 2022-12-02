clear
clc
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library')

%% Step 1: Data
Ym = readmatrix('VS_data','sheet','데이터','range','B6:B152');
X1m = ones(147,1);
X2m = readmatrix('VS_data','sheet','데이터','range','C6:C152');    % 한국 국채 10년물
X3m = readmatrix('VS_data','sheet','데이터','range','D6:D152');    % 국내공급물가지수
X4m = readmatrix('VS_data','sheet','데이터','range','E6:E152');    % 순상품 교역지수
X5m = readmatrix('VS_data','sheet','데이터','range','F6:F152');    % 미국 국채 10년
X6m = readmatrix('VS_data','sheet','데이터','range','H6:H152');    % 종합주가지수
X7m = readmatrix('VS_data','sheet','데이터','range','I6:I152');    % 통화량
X8m = readmatrix('VS_data','sheet','데이터','range','J6:J152');    % 심리지수
Xm = [X1m, X2m, X3m, X4m, X5m, X6m, X7m, X8m];


%% Step 1: Prior
Y = Ym;
X = Xm;
k = cols(X);

stdX = standdc(X);
% sig2 prior
v_0 = 3;
d_0 = 3;

% p prior
c0 = 5;
d0 = 5;

% b0, b1 prior  가장 중요한 prior -> b0와 b1의 크기에 따라 중요한 변수와 중요하지 않은 변수를 선택하는 기준
v_00 = 10;             % prior의 parameter를 변화시키는 것도 bayesian 분석에서 새로운 모형을 적용시키는 것
d_00 = (0.01^2)*10;    % 최적의 prior를 찾는 것도 중요한 연구 중 하나

v_01 = 10;
d_01 = 1*10;

%% Step 2: Sampling
n0 = 1000;    % burn-in
n1 = 1000;    % MCMC size
[betam, sig2m, Gam, postmom_beta, postmom_sig2, postmom_Gam] = MCMC_SSVS(Y, X, n0, n1, v_0, d_0, v_00, d_00, v_01, d_01, c0, d0);

%% Step 3: Results
Gam_hat = postmom_Gam(:, 2);
figure
x = 1:rows(Gam_hat);
scatter(x, Gam_hat, 50, 'b', 'filled')
ylim([-0.05, 1.05])
xlim([0, length(x)+1])
xlabel('Variable')
ylabel('Prob. of inclusion')
grid on
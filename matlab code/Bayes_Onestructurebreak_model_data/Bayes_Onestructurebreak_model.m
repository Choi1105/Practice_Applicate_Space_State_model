%% One Structure Break

clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library')

%% Step 1: Data
Ym = readmatrix('KOSPI','sheet','데이터','range','B2:B110');              % KOSPI 지수 변화율
X1m = ones(109,1);
X2m = readmatrix('KOSPI','sheet','데이터','range','C2:C110');             % PMI 변화율
X3m = readmatrix('KOSPI','sheet','데이터','range','D2:D110'); 


%% Step2: Estimation
% Data
Y = Ym;             % T by 1
X = [X1m, X2m, X3m];     % T by 2

% Information
[T, k] = size(X);

% beta prior
beta_0 = zeros(k,1);
B_0 = 1*eye(k);

% sig2 prior
% mean = d_0/v_0
v_0 = 3;
d_0 = 3;

% uniform for P(probability of break point), T*a<tau<T*b
a_0 = round(0.2*T);
b_0 = round(0.8*T);

%% Step 3: Sampling
n0 = 1000;
n1 = 1000;
n = n0 + n1;
[beta1m, beta2m, sig21m, sig22m, taum, tau_tsm] = Linear_SB(Y, X, beta_0, B_0, v_0, d_0, a_0, b_0, n0, n1);


%% Step 4: Results
% Break points
time = 1:T;
yyaxis left
plot(time, tau_tsm, 'b-', 'LineWidth', 2);
ylabel('Posterior Mass');
ax = gca;
ax.YColor = 'b';

yyaxis right
plot(time, Y, 'k--', 'LineWidth',2);
ax.YColor = 'k';
legend( 'Prob. of breakpoint (LHS)', 'Data (RHS)', 'Location', 'Northwest')

% parameters
Plot_Prior_Post(beta1m(:,1), beta2m(:,1), '$\beta_1$', 'before the break', 'after the braek')
Plot_Prior_Post(beta1m(:,2), beta2m(:,2), '$\beta_2$', 'before the break', 'after the braek')
Plot_Prior_Post(beta1m(:,3), beta2m(:,3), '$\beta_3$', 'before the break', 'after the braek')
Plot_Prior_Post(sqrt(sig21m), sqrt(sig22m), '$\sigma$', 'before the break', 'after the braek')

MHm = [beta1m beta2m sig21m sig22m];
alpha = 0.05;
maxac = 200;
is_plot = 0;
postmom = MHout(MHm, alpha, maxac, is_plot);

% Display
disp('====================================================================');
disp(['    Mean''    S.E.' '      2.5%' '      97.5%' '      ineff' '      p_val']);
disp('====================================================================');
disp([postmom(:,2) postmom(:,3) postmom(:,4) postmom(:,6) postmom(:,7) postmom(:,8)]);
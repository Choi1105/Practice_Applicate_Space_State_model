%% One Structure Break

clear;
clc;
addpath('C:\Users\PC\OneDrive - 전북대학교\문서\MATLAB\M_library')

%% Step 1: DGP
T = 200;               % Sample size
tau = 100;             % 구조변화 시점

% True parameters
b11 = 1;               % 구조변화 이전의 beta 1
b21 = 1;               % 구조변화 이전의 beta 2
b12 = 5;               % 구조변화 이후의 beta 1
b22 = 5;               % 구조변화 이후의 beta 2
sig21 = 0.1;           % 구조변화 이전의 sigma2
sig22 = 0.5;           % 구조변화 이후의 sigma2
tru_para = [b11;b21;b12;b22;sig21;sig22;tau];

Ym = zeros(T,1);       % Y값을 저장할 공간
X1m = ones(T,1);       % T개의 1
X2m = 5*rand(T,1);     % T개의 랜덤값

for t = 1:T            % t이전의 Y1과 Y2

    if t < tau
        Ym(t,1) = b11*X1m(t,1) + b21*X2m(t,1) + randn(1,1)*sqrt(sig21);
    else
        Ym(t,1) = b12*X1m(t,1) + b22*X2m(t,1) + randn(1,1)*sqrt(sig22);
    end

end

%% Step2: Estimation
% Data
Y = Ym;             % T by 1
X = [X1m, X2m];     % T by 2

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
Plot_Prior_Post(sqrt(sig21m), sqrt(sig22m), '$\sigma$', 'before the break', 'after the braek')
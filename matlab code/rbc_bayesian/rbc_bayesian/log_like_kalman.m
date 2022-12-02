%log_like = log_like_kalman(Y,F,Q,H,R)
%Evaluates the the log likelihood function at the sample Y using the Kalman filter. The model is:
% x_t+1 = F x_t + v_t+1
% y_t = H'x_t + w_t 
%var(v_t) = Q
%var(w_t) = R
% x_t is nx by 1
% F is nx by nx
% v_t is nx by 1
% y_t is ny by 1
% H' is ny by nx
% w_t is ny by 1
% Q is nx by nx
% R is ny by ny
%Y is a T by ny matrix containing T observations of y_t
% See Hamilton, page 372
%(c) S. Schmitt_Grohe and M. Uribe, December 11, 2007

function [loglike,lf,s,sP] = log_like_kalman(Y,F,Q,H,R)
% F=hx; Q=nvarshock; H=GX'; R=nvarme;
%Number of Observations
% T = size(Y,1);

%Rows of transition matrix
r = size(F,1);
% n = size(Y,2);

%Initial forecast
% XI_10 = zeros(r,1);

%MSE of initial forecast
%P_10 = zeros(r,r);
%P_10(:) = (eye(r^2) - kron(F,F))\Q(:);
%Doubling algorithm
F_old=F;
Q_old=Q;
P_10_old=eye(size(F));
diferenz=0.1;
while diferenz>1e-25;
    P_10 =F_old*P_10_old*F_old' + Q_old;

    diferenz = max(max(abs(P_10-P_10_old)));
    Q_old=F_old*Q_old*F_old' + Q_old;
    F_old = F_old * F_old;
    P_10_old=P_10;
end    %while diferenz

% XIold = XI_10;
% Pold = P_10;
XIold = zeros(r,1);
Pold = eye(r)*1e-03;


[s,P,eta,f,loglike,lf,sP]=KF_constant(H',F,Q,XIold',Pold,Y',R);

%--------------------------------------------------------------------------
function [s,P,eta,f,llf,lf,sP]=KF_constant(H,F,Q,s00,P00,y,R)
% This function computes the Kalman filter of a simple state space model
% It modifies Prof. Landon-Lane code to include the presence of a 
% drift in the measurement eq.
% y(t)= A + H*s(t) + u_t; Euu' = R,
% s(t)= F*s(t-1) + Q*u(t)
% Output
% s: filtered states (m x T)
% sP: smoothed states (m x T)
% P: P t|t    (m x m x T) 
% eta: forecast error (T x n)
% f: forecast error variance (n x n x T)

[n,T] = size(y); %Note, data is inputed horizontally 
[n,m] = size(H); %Note, m is the number of states 
eta = zeros(T,n); 
%s  = zeros(m,1,T);%Matrix where the updates states are stored. 
 s  = zeros(m,T);
 sp = zeros(m,T);%Matrix where the lagged filtered states are stored. 
P  = zeros(m,m,T); %Matrix where the updated var-cov of states is stored
Pp = zeros(m,m,T); %Matrix where the lagged var-cov of states is stored
lf = zeros(T,1);
l2p = log(2*pi);

%prediction (first obsrvations)
%sp=s00; %Vector of zeros(n,1) initially
 sp(:,1)=s00;
%Pp=P00; 
 Pp(:,:,1)=P00;
%eta(1,:)=(y(:,1)-H*sp)'; %Obtaining first forecast errors 
 eta(1,:)=(y(:,1)-H*sp(:,1))';
%f(:,:,1)=H*Pp*H' + R; 
 f(:,:,1)=H*Pp(:,:,1)*H' + R;    %Note, forecast error variance stored in the first page

%updating (first obsrvations)
%finv= f(:,:,1)\eye(n); 
 finv= pinv(f(:,:,1));
%K=Pp*H'*finv; %Note: NOT the Kalman gain matrix: which is F*Pp*H'*finv
 K=Pp(:,:,1)*H'*finv;
%s(:,1,1)=sp+K*eta(1,:)'; %Note, first filtered states
 s(:,1)=sp(:,1)+K*eta(1,:)';
%P(:,:,1)=(eye(m)-K*H)*Pp; %Note, first filtered var-cov
 P(:,:,1)=(eye(m)-K*H)*Pp(:,:,1);

lf(1)=-(n/2)*l2p+0.5*log(det(finv))-0.5*eta(1,:)*finv*eta(1,:)';

for i=2:T
    %i
    %sp=F*s(:,1,i-1);         %computes St|t-1
     sp(:,i)=F*s(:,i-1);
    %Pp=F*P(:,:,i-1)*F'+Q*Q'; %computes Pt|t-1
     Pp(:,:,i)=F*P(:,:,i-1)*F'+Q;
    %eta(i,:)=(y(:,i)-H*sp)';
     eta(i,:)=(y(:,i)-H*sp(:,i))';
    %f(:,:,i)=H*Pp*H' + R;
     f(:,:,i)=H*Pp(:,:,i)*H' + R;     
    %finv=f(:,:,i)\eye(n);
     finv=pinv(f(:,:,i));
    %K=Pp*H'*finv;
     K=Pp(:,:,i)*H'*finv;
    %s(:,1,i)=sp+K*eta(i,:)';
     s(:,i)=sp(:,i)+K*eta(i,:)';
    %P(:,:,i)=(eye(m)-K*H)*Pp;
     P(:,:,i)=(eye(m)-K*H)*Pp(:,:,i);
    lf(i)=-(n/2)*l2p+0.5*log(det(finv))-0.5*eta(i,:)*finv*eta(i,:)';
end

% compute log likelihood function
llf=sum(lf);

% compute smooth states (uses Hamilton formulae)
stT = s(:,T); % Takes last state updated and takes it as first smoothed
sP(:,T) = stT;
for i = T-1:-1:1 
    J = P(:,:,i)*F'*pinv(Pp(:,:,i+1));
    sP(:,i) = s(:,i) + J*(sP(:,i+1) - sp(:,i+1));
end

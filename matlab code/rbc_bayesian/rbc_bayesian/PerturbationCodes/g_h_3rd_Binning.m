function [gxxxOut,hxxxOut,gssxOut,hssxOut,gsss,hsss] = g_h_3rd_Binning(fx,fxp,fy,fyp,...   % First-order terms
    fxx,fxxp,fxy,fxyp,fxpx,fxpxp,fxpy,fxpyp,...                                 % Second-order terms
    fyx,fyxp,fyy,fyyp,fypx,fypxp,fypy,fypyp,...
    fxxx,fxxxp,fxxy,fxxyp,fxxpxp,fxxpy,fxxpyp,...                                % Third-order terms
    fxyy,fxyyp,fxypyp,fxpxpxp,fxpxpy,fxpxpyp,...
    fxpyy,fxpyyp,fxpypyp,...
    fyyy,fyyyp,fyypyp,fypypyp,...
    gx,gxxSGU,gss,hx,hxxSGU,hss,eta,skew)
% This program computes the matrices gxxx,hxxx,gssx,hssx,gsss,and hsss
% that for a third-order approximation of the DSGE model. That is, if
% E_t[f(yp,y,xp,x)=0, then the solution is of the form
% xp = h(x,sigma) + sigma * eta * ep
% y = g(x,sigma).
%
% The output:
%   - gxxx(beta1,alfa1,alfa2,alfa3) where beta1 = 1,2,...,ny and
%     alfa1,alfa2,alfa3 = 1,2,...nx
%   - hxxx(gama,alfa1,alfa2,alfa3) where gama1 = 1,2,...,nx and
%     alfa1,alfa2,alfa3 = 1,2,...nx
%   - gssx a matrix of order ny * nx
%   - hssx a matrix of order nx * nx
%   - gsss a matrix of order ny * 1
%   - hsss a matrix of order nx * 1
% The input:
%   - skew is a matrix of dimension nx * nx^2.
%     skew = E_t[u_t kron u_t' kron u_t'] where u_t = sig*eta*eta_t;
%   - all other inputs are standard ...
%
% Copyright Andrew Binning 2011.
% Remark: These codes are a modification of third_order.m by Binning, to
% allow for the same input arguments as used in the codes by Andreasen in
% g_h_3rd_v1 and g_h_3rd_v2

sigma = eta*eta';
nx = size(hx,1);
ny = size(gx,1);
n = nx + ny;

% By MMA: Computing D,H,T
[D,H,T] = getDHT(n,nx,ny,fx,fxp,fy,fyp,...
    fxx,fxxp,fxy,fxyp,fxpx,fxpxp,fxpy,fxpyp,...
    fyx,fyxp,fyy,fyyp,fypx,fypxp,fypy,fypyp,...
    fxxx,fxxxp,fxxy,fxxyp,fxxpxp,fxxpy,fxxpyp,...
    fxyy,fxyyp,fxypyp,fxpxpxp,fxpxpy,fxpxpyp,...
    fxpyy,fxpyyp,fxpypyp,...
    fyyy,fyyyp,fyypyp,fypypyp);

% By MMA: Constructing hxxa and gxxa from the representation of hxx and gxx
% used by SGU.
hxxa = reshape(hxxSGU,nx,nx^2);
gxxa = reshape(gxxSGU,ny,nx^2);

% By MMA: Constructing hxx and gxx as used by Binning
hxx  = zeros(nx^2,nx);
gxx  = zeros(ny*nx,nx);
for i=1:nx
    hxx((i-1)*nx+1:i*nx,:) = squeeze(hxxSGU(i,:,:));
end
for i=1:ny
    gxx((i-1)*nx+1:i*nx,:) = squeeze(gxxSGU(i,:,:));
end

%--------------------------------------------------------------------------
% Solve the non-stochastic block of the model
%--------------------------------------------------------------------------

MX = [hx;gx*hx;eye(nx);gx];
MXX = [hxx;kron3(eye(ny),hx')*gxx*hx + kron3(gx,eye(nx))*hxx;zeros(nx^2,nx);gxx];
MXXA = [hxxa;gxxa*kron3(hx,hx) + gx*hxxa;zeros(nx,nx^2);gxxa];

d1 = D(:,1:nx);
d2 = D(:,nx+1:n);
d4 = D(:,n+nx+1:2*n);

a1 = kron3(eye(ny),kron3(hx',eye(nx)))*kron3(gxx,eye(nx))*hxx;

temp_hx = nan(nx^2,nx^2);

for ii = 1:nx
    temp_hx(nx*(ii-1)+(1:nx),:) = kron3(eye(nx),hx(ii,:));
end

a2 = kron3(eye(ny),temp_hx')*kron3(gxx,eye(nx))*hxx;
a3 = kron3(eye(ny),hxxa')*gxx*hx;

eq1 = kron3(eye(n),kron3(MX',MX'))*T*MX;
eq2 = kron3(eye(n),MXXA')*H*MX;
eq3 = kron3(eye(n),kron3(MX',eye(nx)))*kron3(H,eye(nx))*MXX;

temp_MX = nan(2*n*nx,nx^2);

for ii = 1:2*n
    temp_MX(nx*(ii-1)+(1:nx),:) = kron3(eye(nx),MX(ii,:));
end

eq4 = kron3(eye(n),temp_MX')*kron3(H,eye(nx))*MXX;

eq5 = kron3(d2,eye(nx^2))*(a1+a2+a3);

A1 = eq1 + eq2 + eq3 + eq4 + eq5;

B1 = kron3(d1,eye(nx^2));
B2 = kron3(d2,eye(nx^2))*kron3(eye(ny),kron3(hx',hx'));
B3 = kron3(d2,eye(nx^2))*kron3(gx,eye(nx^2));
B4 = kron3(d4,eye(nx^2));

lhs = [kron3(eye(nx),B1+B3),kron3(hx',B2)+kron3(eye(nx),B4)];
rhs = reshape(-A1,[],1);

result1 = lhs\rhs;
hxxx    = reshape(result1(1:nx^4),nx^3,nx);
gxxx    = reshape(result1(nx^4+1:end),ny*nx^2,nx);

% By MMA: representing the solution in multi-arrays
hxxxOut = zeros(nx,nx,nx,nx);
p       = 0;
for ii = 1:nx
    for jj = 1:nx
        p = p + 1;
        hxxxOut(ii,:,:,jj) = hxxx((1:nx)+nx*(p-1),:);
    end
end
gxxxOut = zeros(ny,nx,nx,nx);
p       = 0;
for ii = 1:ny
    for jj = 1:nx
        p = p + 1;
        gxxxOut(ii,:,:,jj) = gxxx((1:nx)+nx*(p-1),:);
    end
end


%--------------------------------------------------------------------------
% Solve for the risk adjusted to the slope terms
%--------------------------------------------------------------------------

N = [eye(nx);gx;zeros(n,nx)];
NXXA = [zeros(nx,nx^2);gxxa*kron3(hx,eye(nx));zeros(n,nx^2)];

BB = [hss;tracem(kron3(eye(ny),sigma)*gxx) + gx*hss + gss; zeros(nx,1); gss];

req1 = tracem(kron3(eye(n),kron3(MX',N'))*T*N*sigma);
req2 = tracem(kron3(eye(n),NXXA')*H*N*sigma);
req3 = kron3(eye(n),MX')*H*BB;

CC = tracem(kron3(eye(nx*ny),sigma)*kron3(eye(ny),kron3(hx',eye(nx)))*gxxx)+kron3(eye(ny),hx')*gxx*hss;

req6 = kron3(d2,eye(nx))*CC;

A2 = req1 + 2*req2 + req3 + req6;


BB = [kron3(d1,eye(nx)) + kron3(d2,eye(nx))*kron3(gx,eye(nx)),kron3(d2,hx') + kron3(d4,eye(nx))];

% Note that kron3(d2,hx') = kron3(d2,eye(nx))*kron3(eye(ny),hx')

result2 = -BB\A2;

hssx = result2(1:nx^2,:);
gssx = result2(nx^2+1:end);

% By MMA: representing the solution in a matrix
hssxOut = reshape(hssx,nx,nx)';
gssxOut = reshape(gssx,nx,ny)';

%--------------------------------------------------------------------------
% Solve for the effects of skewness
%--------------------------------------------------------------------------

if all(skew == 0)
    hsss = zeros(nx,1);
    gsss = zeros(ny,1);
else
    NSSA = [zeros(nx,nx^2);gxxa;zeros(n,nx^2)];
    
    xeq1 = tracem(kron3(eye(n),kron3(N',N'))*T*N*skew);
    xeq2 = tracem(kron3(eye(n),NSSA')*H*N*skew);
    xeq3 = d2*tracem(kron3(eye(ny),skew)*gxxx);
    
    A3 = xeq1 + 3*xeq2 + xeq3;
    
    GG = [d1 + d2*gx, d2 + d4];
    
    result3 = -GG\A3;
    
    hsss = result3(1:nx,:);
    gsss = result3(nx+1:end,:);
end

end

function [D,H,T] = getDHT(n,nx,ny,fx,fxp,fy,fyp,...
    fxx,fxxp,fxy,fxyp,fxpx,fxpxp,fxpy,fxpyp,...
    fyx,fyxp,fyy,fyyp,fypx,fypxp,fypy,fypyp,...
    fxxx,fxxxp,fxxy,fxxyp,fxxpxp,fxxpy,fxxpyp,...
    fxyy,fxyyp,fxypyp,fxpxpxp,fxpxpy,fxpxpyp,...
    fxpyy,fxpyyp,fxpypyp,...
    fyyy,fyyyp,fyypyp,fypypyp)
% This function is written by MMA, by adopting some of Binnings codes.

%--------------------------------------------------------------------------
% first derivatives
%--------------------------------------------------------------------------
D = [fxp,fyp,fx,fy];

%--------------------------------------------------------------------------
% Second derivatives
%--------------------------------------------------------------------------
Htmp  = zeros(n,2*n,2*n);
xpr = 1:nx;
ypr = nx+1:n;
xr  = n+1:n+nx;
yr  = n+nx+1:2*n;

Htmp(:,xpr,xpr) = fxpxp;
Htmp(:,xpr,ypr) = fxpyp;
Htmp(:,xpr,yr)  = fxpy;
Htmp(:,xpr,xr)  = fxpx;
Htmp(:,ypr,xpr) = fypxp;
Htmp(:,ypr,ypr) = fypyp;
Htmp(:,ypr,xr)  = fypx;
Htmp(:,ypr,yr)  = fypy;
Htmp(:,xr,xpr)  = fxxp;
Htmp(:,xr,ypr)  = fxyp;
Htmp(:,xr,xr)   = fxx;
Htmp(:,xr,yr)   = fxy;
Htmp(:,yr,xpr)  = fyxp;
Htmp(:,yr,ypr)  = fyyp;
Htmp(:,yr,xr)   = fyx;
Htmp(:,yr,yr)   = fyy;

H               = zeros(2*n^2,2*n);
for ii = 1:n
    for jj = 1:2*n
        H(jj+(ii-1)*2*n,:) = Htmp(ii,:,jj);
    end
end
%--------------------------------------------------------------------------
% Third derivatives
%--------------------------------------------------------------------------
% We start by computing the remaining matrices for the third order derivatives.
% 1) Derivatives related to fxx
%fxxx  : is computed
%fxxxp : is computed
%fxxy  : is computed
%fxxyp : is computed

% 2) Derivatives related to fxxp
fxxpx  = permute(fxxxp,[1,2,4,3]);
%fxxpxp  : is computed
%fxxpy   : is computed
%fxxpyp  : is computed

% 3) Derivatives related to fxy
fxyx  = permute(fxxy,[1,2,4,3]);
fxyxp = permute(fxxpy,[1,2,4,3]);
%fxyy    : is computed
%fxyyp   : is computed

% 4) Derivatives related to fxyp
fxypx  = permute(fxxyp,[1,2,4,3]);
fxypxp = permute(fxxpyp,[1,2,4,3]);
fxypy  = permute(fxyyp,[1,2,4,3]);
%fxypyp   : is computed

% 5) Derivatives related to fxpx
fxpxx  = permute(fxxpx,[1,3,2,4]);
fxpxxp = permute(fxxpxp,[1,3,2,4]);
fxpxy  = permute(fxxpy,[1,3,2,4]);
fxpxyp = permute(fxxpyp,[1,3,2,4]);

% 6) Derivatives related to fxpxp
fxpxpx = permute(fxxpxp,[1,3,4,2]);
%fxpxpxp   : is computed
%fxpxpy    : is computed
%fxpxpyp   : is computed

% 7) Derivatives related to fxpy
fxpyx  = permute(fxxpy,[1,3,4,2]);
fxpyxp = permute(fxpxpy,[1,2,4,3]);
%fxpyy   : is computed
%fxpyyp  : is computed

% 8) Derivatives related to fxpyp
fxpypx  = permute(fxypxp,[1,4,3,2]);
fxpypxp = permute(fxpxpyp,[1,2,4,3]);
fxpypy  = permute(fxpyyp,[1,2,4,3]);
%fxpyyp : is computed

% 9) Derivatives related to fyx
fyxx  = permute(fxxy,[1,4,2,3]);
fyxxp = permute(fxpxy,[1,4,3,2]);
fyxy  = permute(fxyy,[1,3,2,4]);
fyxyp = permute(fxyyp,[1,3,2,4]);

% 10) Derivatives related to fyxp
fyxpx  = permute(fxxpy,[1,4,3,2]);
fyxpxp = permute(fxpxpy,[1,4,2,3]);
fyxpy  = permute(fxpyy,[1,3,2,4]);
fyxpyp = permute(fxpyyp,[1,3,2,4]);

% 11) Derivatives related to fyy
fyyx  = permute(fxyy,[1,3,4,2]);
fyyxp = permute(fxpyy,[1,3,4,2]);
%fyyy  : is computed
%fyyyp : is computed

% 12) Derivatives related to fyyp
fyypx  = permute(fxyyp,[1,3,4,2]);
fyypxp = permute(fxpyyp,[1,3,4,2]);
fyypy  = permute(fyyyp,[1,2,4,3]);
%fyypyp  : is computed

% 13) Derivatives related to fypx
fypxx  = permute(fxxyp,[1,4,2,3]);
fypxxp = permute(fxxpyp,[1,4,2,3]);
fypxy  = permute(fxyyp,[1,4,2,3]);
fypxyp = permute(fxypyp,[1,3,2,4]);

% 14) Derivatives related to fypxp
fypxpx  = permute(fxxpyp,[1,4,3,2]);
fypxpxp = permute(fxpxpyp,[1,4,2,3]);
fypxpy  = permute(fxpyyp,[1,4,2,3]);
fypxpyp = permute(fxpypyp,[1,3,2,4]);

% 15) Derivatives related to fypy
fypyx  = permute(fxyyp,[1,4,3,2]);
fypyxp = permute(fxpyyp,[1,4,3,2]);
fypyy  = permute(fyyyp,[1,4,2,3]);
fypyyp = permute(fyypyp,[1,4,2,3]);

% 16) Derivatives related to fypyp
fypypx = permute(fxypyp,[1,3,4,2]);
fypypxp = permute(fxpypyp,[1,3,4,2]);
fypypy = permute(fyypyp,[1,3,4,2]);
%fypypyp  : is computed

Ttmp = zeros(n,2*n,2*n,2*n);
Ttmp(:,xpr,xpr,xpr) = fxpxpxp;
Ttmp(:,xpr,ypr,xpr) = fxpypxp;
Ttmp(:,xpr,yr,xpr)  = fxpyxp;
Ttmp(:,xpr,xr,xpr)  = fxpxxp;
Ttmp(:,ypr,xpr,xpr) = fypxpxp;
Ttmp(:,ypr,ypr,xpr) = fypypxp;
Ttmp(:,ypr,xr,xpr)  = fypxxp;
Ttmp(:,ypr,yr,xpr)  = fypyxp;
Ttmp(:,xr,xpr,xpr)  = fxxpxp;
Ttmp(:,xr,ypr,xpr)  = fxypxp;
Ttmp(:,xr,xr,xpr)   = fxxxp;
Ttmp(:,xr,yr,xpr)   = fxyxp;
Ttmp(:,yr,xpr,xpr)  = fyxpxp;
Ttmp(:,yr,ypr,xpr)  = fyypxp;
Ttmp(:,yr,xr,xpr)   = fyxxp;
Ttmp(:,yr,yr,xpr)   = fyyxp;

Ttmp(:,xpr,xpr,ypr) = fxpxpyp;
Ttmp(:,xpr,ypr,ypr) = fxpypyp;
Ttmp(:,xpr,yr,ypr)  = fxpyyp;
Ttmp(:,xpr,xr,ypr)  = fxpxyp;
Ttmp(:,ypr,xpr,ypr) = fypxpyp;
Ttmp(:,ypr,ypr,ypr) = fypypyp;
Ttmp(:,ypr,xr,ypr)  = fypxyp;
Ttmp(:,ypr,yr,ypr)  = fypyyp;
Ttmp(:,xr,xpr,ypr)  = fxxpyp;
Ttmp(:,xr,ypr,ypr)  = fxypyp;
Ttmp(:,xr,xr,ypr)   = fxxyp;
Ttmp(:,xr,yr,ypr)   = fxyyp;
Ttmp(:,yr,xpr,ypr)  = fyxpyp;
Ttmp(:,yr,ypr,ypr)  = fyypyp;
Ttmp(:,yr,xr,ypr)   = fyxyp;
Ttmp(:,yr,yr,ypr)   = fyyyp;

Ttmp(:,xpr,xpr,xr)  = fxpxpx;
Ttmp(:,xpr,ypr,xr)  = fxpypx;
Ttmp(:,xpr,yr,xr)   = fxpyx;
Ttmp(:,xpr,xr,xr)   = fxpxx;
Ttmp(:,ypr,xpr,xr)  = fypxpx;
Ttmp(:,ypr,ypr,xr)  = fypypx;
Ttmp(:,ypr,xr,xr)   = fypxx;
Ttmp(:,ypr,yr,xr)   = fypyx;
Ttmp(:,xr,xpr,xr)   = fxxpx;
Ttmp(:,xr,ypr,xr)   = fxypx;
Ttmp(:,xr,xr,xr)    = fxxx;
Ttmp(:,xr,yr,xr)    = fxyx;
Ttmp(:,yr,xpr,xr)   = fyxpx;
Ttmp(:,yr,ypr,xr)   = fyypx;
Ttmp(:,yr,xr,xr)    = fyxx;
Ttmp(:,yr,yr,xr)    = fyyx;

Ttmp(:,xpr,xpr,yr)  = fxpxpy;
Ttmp(:,xpr,ypr,yr)  = fxpypy;
Ttmp(:,xpr,yr,yr)   = fxpyy;
Ttmp(:,xpr,xr,yr)   = fxpxy;
Ttmp(:,ypr,xpr,yr)  = fypxpy;
Ttmp(:,ypr,ypr,yr)  = fypypy;
Ttmp(:,ypr,xr,yr)   = fypxy;
Ttmp(:,ypr,yr,yr)   = fypyy;
Ttmp(:,xr,xpr,yr)   = fxxpy;
Ttmp(:,xr,ypr,yr)   = fxypy;
Ttmp(:,xr,xr,yr)    = fxxy;
Ttmp(:,xr,yr,yr)    = fxyy;
Ttmp(:,yr,xpr,yr)   = fyxpy;
Ttmp(:,yr,ypr,yr)   = fyypy;
Ttmp(:,yr,xr,yr)    = fyxy;
Ttmp(:,yr,yr,yr)    = fyyy;

T                   = zeros(4*n^3,2*n);
for ii = 1:n
    for jj = 1:2*n
        for kk = 1:2*n
            T(kk + (jj-1)*2*n +(ii-1)*4*n^2,:) = Ttmp(ii,:,jj,kk);
        end
    end
end

end


function y = tracem(x)
n = size(x,2);
m = size(x,1)/n;
y = zeros(m,1);
for i=1:m
    y(i,1)=trace(x((n*(i-1)+1):i*n,1:n));
end

end

function Z=kron3(X,Y)
%KRON   Kronecker tensor product. 
Z=reshape(permute(reshape(Y(:)*X(:).',[size(Y),size(X)]),[1,3,2,4]),[size(Y).*size(X)]);
end

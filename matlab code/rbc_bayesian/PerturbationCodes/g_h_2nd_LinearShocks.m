% Copyright Andrew Binning (2013)

function [gxxOut,hxxOut,hss,gss] = g_h_2nd_LinearShocks(D,H,gx,hx,ssigma,nx,ny,n,mx,my,myx);
% EQUATIONS:
% The first m = mx + my equations must be the endogenous equations.
%
% STATE VARIABLES
% The first mx state variables must be endogenous and only predetermined (e.g. capital).
% The next myx variables must be both predetermined and non-predetermined 
% (e.g. lagged consumption). 
% The last nx - mx - myx variables must be exogenous linear processes.  
%
% CONTROL VARIABLES
% The first myx control variables must be both predetermined and
% non-predetermined (e.g. consumption), and with the same ordering as in
% the state variables.
% The next my - myx control variables must be endogenous and only non-predetermined.  
% The last ny - my control variables will be exogenous variables that occur in period t + 1 
% (e.g. expected consumption growth). 
%
%
% INFORMATION
% If one lets mx = nx, myx = 0 and my = ny, then one runs the code without taking into account the linearity 
% of the shock processes and the doubling up of lagged nonpredetermined variables. 
%
% News:
% 9/8 2014: updated to allow lagged control variables (e.g. lagged
% consumption) leaded one period to also affect the model. This implies
% that we may have c_ba1p in the codes instead of c_cu. 

f1 = D(:,1:nx);
f2 = D(:,nx+1:n);
f4 = D(:,n+nx+1:2*n);

MX = [hx;gx*hx;eye(nx);gx];

%syldd = H(1:mx+my,:)*kron3(MX,MX);

syldd = nan(mx+my,nx^2);

for ii = 1:mx+my
    syldd_tmp = KronTimesVector(MX',MX',H(ii,:)')';
    syldd(ii,:) = syldd_tmp;
end

%sylaa = [f1+f2*gx,f4];
sylaa = [f1(1:mx+my,1:mx)+f2(1:mx+my,:)*gx(:,1:mx),f1(1:mx+my,mx+1:mx+myx) + f4(1:mx+my,1:myx) + f2(1:mx+my,:)*gx(:,mx+1:mx+myx), f4(1:mx+my,myx+1:my)];

%sylbb = [zeros(n,nx),f2];
sylbb = [zeros(mx+my,mx),f2(1:mx+my,1:myx),f2(1:mx+my,myx+1:my)];

temp = sylaa\sylbb;
AA = -sylaa\syldd;

clear sylaa slybb syldd

sylcc = hx;

[U,K] = schur(temp);
[V,F] = schur(sylcc);

Dbar = U'*AA*kron(V,V);

clear AA

d = Dbar(:);

clear Dbar

r = 1;

k = 2;

y = solv1(r,F',K,d,k);

clear F K d

xx = U*reshape(y,mx+my,nx^2)*kron3(V',V');

clear U y V 

hxx = [xx(1:mx+myx,:);zeros(nx-(mx+myx),nx^2)];
gxx = [xx(mx+1:mx+my,:);zeros(ny-my,nx^2)];
% By MMA:
hxxOut = reshape(hxx,nx,nx,nx);
gxxOut = reshape(gxx,ny,nx,nx);

clear xx

%lhs = [f1+f2*gx, f2+f4];

lhs = [f1(1:mx+my,1:mx) + f2(1:mx+my,:)*gx(:,1:mx), f2(1:mx+my,1:myx) + f4(1:mx+my,1:myx) + f1(1:mx+my,mx+1:mx+myx) + f2(1:mx+my,:)*gx(:,mx+1:mx+myx),...
    f2(1:mx+my,myx+1:my) + f4(1:mx+my,myx+1:my)];

NX = [eye(nx);gx;zeros(n,nx)];


%tmp1 = gxx*kron3(eye(nx),ssigma);

tmp1 = nan(ny,nx^2);

for ii = 1:ny
    tmp1_tmp = KronTimesVector(eye(nx),ssigma',gxx(ii,:)')';
    tmp1(ii,:) = tmp1_tmp;
end

clear tmp1_tmp

%tmp2 = H(1:mx+my,:)*kron(NX,NX*ssigma);

tmp2 = nan(ny,nx^2);

for ii = 1:mx+my
    tmp2_tmp = KronTimesVector(NX',ssigma'*NX',H(ii,:)')';
    tmp2(ii,:) = tmp2_tmp;
end

clear tmp2_tmp

tmp_tr1 = zeros(ny,1);
tmp_tr2 = zeros(mx+my,1);

for ii = 1:nx
    tmp_tr1 = tmp1(:,ii+(ii-1)*nx) + tmp_tr1;
    tmp_tr2 = tmp2(:,ii+(ii-1)*nx) + tmp_tr2;
end

clear tmp1 tmp2

rhs = f2(1:mx+my,:)*tmp_tr1 + tmp_tr2;

clear tmp_tr1 tmp_tr2

ss = -lhs\rhs;


hss = [ss(1:mx+myx);zeros(nx-mx-myx,1)];
gss = [ss(mx+1:mx+my);zeros(ny-my,1)];

clear ss
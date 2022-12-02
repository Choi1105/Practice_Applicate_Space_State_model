% Copyright Andrew Binning 2013

function [gxxxOut,hxxxOut,gssx,hssx,gsss,hsss] = g_h_3rd_LinearShocks(D,H,T,gx,hx,gxxSGU,hxxSGU,gss,hss,sigma,skew,nx,ny,n,mx,my,myx)
% EQUATIONS:
% The first m = mx + my equations must be the endogenous equations.
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
% INFORMATION
% If one lets mx = nx, myx = 0 and my = ny, then one runs the code without taking into account the linearity 
% of the shock processes and the doubling up of lagged nonpredetermined variables. 
%
% News:
% 9/8 2014: updated to allow lagged control variables (e.g. lagged
% consumption) leaded one period to also affect the model. This implies
% that we may have c_ba1p in the codes instead of c_cu.

%--------------------------------------------------------------------------
% Solve the non-stochastic block of the model
%--------------------------------------------------------------------------

% By MMA: Constructing hxxa and gxxa from the representation of hxx and gxx
% used by SGU.
hxx = reshape(hxxSGU,nx,nx^2);
gxx = reshape(gxxSGU,ny,nx^2);

MX = [hx;gx*hx;eye(nx);gx];
MXXA = [hxx;gxx*kron(hx,hx) + gx*hxx;zeros(nx,nx^2);gxx];

f1 = D(:,1:nx);
f2 = D(:,nx+1:n);
f4 = D(:,n+nx+1:2*n);

a1 = gxx*kron3(hx,hxx);

temp_hxx = nan(nx^2,nx^3);

for ii = 1:nx
    
    temp_hxx(:,1+(ii-1)*nx^2:ii*nx^2) = kron3(eye(nx), hxx(:,1+(ii-1)*nx:ii*nx));
    
end

a2 = gxx*kron3(hx,eye(nx))*temp_hxx;

clear temp_hxx

a3 = gxx*kron3(hxx,hx);

%eq1 = T*kron(MX,kron(MX,MX));

eq1 = nan(mx+my,nx^3);

tmp_MXMX = kron3(MX',MX');

for ii = 1:mx+my
    eq1_tmp = KronTimesVector(tmp_MXMX,MX',T(ii,:)')';
    eq1(ii,:) = eq1_tmp;
end

clear tmp_MXMX

% eq2 = H(1:mx+my,:)*kron3(MX,MXXA);

eq2 = nan(mx+my,nx^3);

for ii = 1:mx+my
    eq2_tmp = KronTimesVector(MX',MXXA',H(ii,:)')';
    eq2(ii,:) = eq2_tmp;
end

clear eq2_tmp

MXX_temp = nan(2*n*nx,nx^3);

for ii = 1:nx
    
    tmp = zeros([2*n, nx, nx, nx]);
    tmp(:,:,1:nx+1:nx^2) = repmat(MXXA(:,1+(ii-1)*nx:ii*nx), [1, 1, nx]);
    tmp = permute(tmp, [1 3 2 4]);
    tmp = reshape(tmp, nx*[2*n, nx]);
    
    MXX_temp(:,1+(ii-1)*nx^2:ii*nx^2) = tmp;
    
    %MXX_temp(:,1+(ii-1)*nx^2:ii*nx^2) = kron(eye(nx),MXXA(:,1+(ii-1)*nx:ii*nx));
    
end

clear tmp

eq3 = H(1:mx+my,:)*kron3(MX,eye(2*n))*MXX_temp;

clear MXX_temp

%eq4 = H(1:mx+my,:)*kron3(MXXA,MX);

eq4 = nan(mx+my,nx^3);

for ii = 1:mx+my
    eq4_tmp = KronTimesVector(MXXA',MX',H(ii,:)')';
    eq4(ii,:) = eq4_tmp;
end

clear eq4_tmp

eq5 = f2(1:mx+my,:)*(a1+a2+a3);

clear a1 a2 a3

syldd = eq1 + eq2 + eq3 + eq4 + eq5;

clear eq1 eq2 eq3 eq4 eq5

sylaa = [f1(1:mx+my,1:mx) + f2(1:mx+my,:)*gx(:,1:mx), f1(1:mx+my,mx+1:mx+myx) + f2(1:mx+my,:)*gx(:,mx+1:mx+myx) + f4(1:mx+my,1:myx), f4(1:mx+my,myx+1:my)];

sylbb = [zeros(mx+my,mx),f2(1:mx+my,1:my)];

temp = sylaa\sylbb;
AA = -sylaa\syldd;

clear sylaa slybb syldd

[U,K] = schur(temp);
[V,F] = schur(hx);

%keyboard

%Dbar1 = U'*AA*kron3(kron3(V,V),V);

Dbar = nan(mx+my,nx^3);

tmp_UAA = AA'*U;
tmp_VV = kron3(V',V');

for ii = 1:mx+my
    Dbar(ii,:) = KronTimesVector(tmp_VV,V',tmp_UAA(:,ii))';
end

clear tmp_UAA tmp_VV

clear AA

d = Dbar(:);

clear Dbar

r = 1;

k = 3;

y = solv1(r,F',K,d,k);

clear F K d

%xxx = U*reshape(y,mx+my,nx^3)*kron(kron(V',V'),V');

tmp_VV = kron3(V,V);

tmp_Uy = reshape(y,mx+my,nx^3)'*U';

clear y U

xxx = nan(mx+my,nx^3);

for ii = 1:mx+my
    xxx(ii,:) = KronTimesVector(tmp_VV,V,tmp_Uy(:,ii))';
end

clear tmp_VV V tmp_Uy

hxxx = [xxx(1:mx+myx,:);zeros(nx-mx-myx,nx^3)];
gxxx = [xxx(mx+1:mx+my,:);zeros(ny-my,nx^3)];

clear xxx

% By MMA:
hxxxOut = reshape(hxxx,nx,nx,nx,nx);
gxxxOut = reshape(gxxx,ny,nx,nx,nx);

clear xxx

%--------------------------------------------------------------------------
% Solve for the risk adjusted to the slope terms
%--------------------------------------------------------------------------

N = [eye(nx);gx;zeros(n,nx)];
NXXA = [zeros(nx,nx^2);gxx*kron(hx,eye(nx));zeros(n,nx^2)];

tmp1 = gxx*kron3(eye(nx),sigma);

tmp_tr1 = zeros(ny,1);

for ii = 1:nx
    tmp_tr1 = tmp1(:,ii+(ii-1)*nx) + tmp_tr1;
end

BB = [hss;tmp_tr1 + gx*hss + gss; zeros(nx,1); gss];

clear tmp_tr1

%req1 = tracem(kron(eye(n),kron(MX',N'))*T*N*sigma);

tmp2 = nan(mx+my,nx^3);

tmp_MXN = kron3(MX',N');
tmp_sigN = sigma'*N';

for ii = 1:mx+my
    tmp2_tmp = KronTimesVector(tmp_MXN,tmp_sigN,T(ii,:)')';
    tmp2(ii,:) = tmp2_tmp;
end

clear tmp_MXN tmp_sigN

%tmp2 = T*kron3(kron3(MX,N),N)*kron3(eye(nx^2),sigma);

%tmp3 = H(1:mx+my,:)*kron3(NXXA,N*sigma);

tmp3 = nan(mx+my,nx^3);

for ii = 1:mx+my
    tmp3_tmp = KronTimesVector(NXXA',sigma'*N',H(ii,:)')';
    tmp3(ii,:) = tmp3_tmp;
end

clear tmp3_tmp

%tmp4 = gxxx*kron3(hx,eye(nx^2))*kron3(eye(nx^2),sigma);

tmp4_tmp1 = kron(sparse(eye(nx^2)),sparse(sigma));

tmp4_tmp2 = kron(hx,sparse(eye(nx^2)));

tmp4 = gxxx*tmp4_tmp2*tmp4_tmp1;

clear tmp4_tmp2 tmp4_tmp1

req1 = nan(mx+my,nx);
req2 = nan(mx+my,nx);
CC = nan(ny,nx);

for ii = 1:nx
    tmp_tr2 = zeros(mx+my,1);
    tmp_tr3 = zeros(mx+my,1);
    tmp_tr4 = zeros(ny,1);
    for jj = 1:nx
        tmp_tr2 = tmp2(:,jj+(jj-1)*nx + (ii-1)*nx^2) + tmp_tr2;
        tmp_tr3 = tmp3(:,jj+(jj-1)*nx + (ii-1)*nx^2) + tmp_tr3;
        tmp_tr4 = tmp4(:,jj+(jj-1)*nx + (ii-1)*nx^2) + tmp_tr4;
    end
    req1(:,ii) = tmp_tr2;
    req2(:,ii) = tmp_tr3;
    CC(:,ii) = tmp_tr4;
end

clear tmp2 tmp3 tmp4 tmp_tr2 tmp_tr3 tmp_tr4

req3 = H(1:mx+my,:)*kron3(MX,BB);

clear BB NXXA MX

CC = CC + gxx*kron3(hx,hss);

req6 = f2(1:mx+my,:)*CC;

clear CC

syldd = req1 + 2*req2 + req3 + req6;

clear req1 req2 req3 req6

sylaa = [f1(1:mx+my,1:mx) + f2(1:mx+my,:)*gx(:,1:mx), f1(1:mx+my,mx+1:mx+myx) + f2(1:mx+my,:)*gx(:,mx+1:mx+myx) + f4(1:mx+my,1:myx), f4(1:mx+my,myx+1:my)];

sylbb = [zeros(mx+my,mx),f2(1:mx+my,1:my)];

temp = sylaa\sylbb;
AA = -sylaa\syldd;

clear sylaa slybb syldd

[U,K] = schur(temp);
[V,F] = schur(hx);

Dbar = U'*AA*V;

clear AA

d = Dbar(:);

clear Dbar

r = 1;

k = 1;

y = solv1(r,F',K,d,k);

clear F K d

xx = U*reshape(y,mx+my,nx)*V';

clear U y V

hssx = [xx(1:mx+myx,:);zeros(nx-mx-myx,nx)];
gssx = [xx(mx+1:mx+my,:);zeros(ny-my,nx)];

clear xx

%--------------------------------------------------------------------------
% Solve for the effects of skewness
%--------------------------------------------------------------------------
if all(skew == 0)
    hsss = zeros(nx,1);
    gsss = zeros(ny,1);
else
    NSSA = [zeros(nx,nx^2);gxx;zeros(n,nx^2)];
    
    %xeq1 = T*kron(kron(N,N),N*skew);
    
    tmp5 = nan(mx+my,nx^4);
    
    %temp_N = kron3(skew'*N',N');
    temp_N = kron(sparse(skew'*N'),N');
    
    for ii = 1:mx+my
        tmp5_tmp = KronTimesVector(temp_N,N',T(ii,:)')';
        tmp5(ii,:) = tmp5_tmp;
    end
    
    clear temp_N
    
    %tmp6 = H*kron3(NSSA,N*skew);
    
    tmp6 = nan(mx+my,nx^4);
    
    for ii = 1:mx+my
        tmp6_tmp = KronTimesVector(NSSA',sparse(skew'*N'),H(ii,:)')';
        tmp6(ii,:) = tmp6_tmp;
    end
    
    tmp7 = gxxx*kron(sparse(eye(nx^2)),skew);
    
    xeq1 = zeros(mx+my,1);
    xeq2 = zeros(mx+my,1);
    tmp7_tmp = zeros(ny,1);
    
    for ii = 1:nx^2
        xeq1 = tmp5(:,ii+(ii-1)*nx^2) + xeq1;
        xeq2 = tmp6(:,ii+(ii-1)*nx^2) + xeq2;
        tmp7_tmp = tmp7(:,ii+(ii-1)*nx^2) + tmp7_tmp;
    end
    
    clear tmp5 tmp6 tmp7
    
    xeq3 = f2(1:mx+my,:)*tmp7_tmp;
    
    clear tmp7_tmp
    
    RRR = xeq1 + 3*xeq2 + xeq3;
    
    clear xeq1 xeq2 xeq3
    
    GGG = [f1(1:mx+my,1:mx) + f2(1:mx+my,:)*gx(:,1:mx),...
        f1(1:mx+my,mx+1:mx+myx) + f2(1:mx+my,:)*gx(:,mx+1:mx+myx) + f2(1:mx+my,1:myx) + f4(1:mx+my,1:myx),...
        f2(1:mx+my,myx+1:my) + f4(1:mx+my,myx+1:my)];
    
    result3 = -GGG\RRR;
    
    clear GGG RRR
    
    hsss = [result3(1:mx+myx);zeros(nx-mx-myx,1)];
    gsss = [result3(mx+1:end,:);zeros(ny-my,1)];
    
end

return


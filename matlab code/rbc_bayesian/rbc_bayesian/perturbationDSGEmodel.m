% By Sanha Noh 2020.03.10

function [gx,hx,gxx,hxx,gss,hss,gxxx,hxxx,gssx,hssx,gsss,hsss,eta,n,ny,nx,ne] = perturbationDSGEmodel(params,order_app,linearShocksOn,logApprox);
         
%% Computing the steady state values
% Number of equations in the model

[nf,nfx,nfxp,nfy,nfyp,eta,errorMes] = numDerivF_1st(params);
if order_app > 1
    [nfxx,nfxxp,nfxy,nfxyp,nfxpx,nfxpxp,nfxpy,nfxpyp,nfyx,nfyxp,nfyy,nfyyp,nfypx,nfypxp,nfypy,nfypyp,errorMes] = numDerivF_2nd(params);
end
if order_app > 2
    [nfxxx,nfxxxp,nfxxy,nfxxyp,nfxxpxp,nfxxpy,nfxxpyp,nfxyy,nfxyyp,nfxypyp,nfxpxpxp,nfxpxpy,nfxpxpyp,nfxpyy,nfxpyyp,nfxpypyp,nfyyy,nfyyyp,nfyypyp,nfypypyp,errorMes] = numDerivF_3rd(params);
end
nx = size(nfx,2);
ny = size(nfy,2);
ne = size(eta,2);
n = ny + nx;

%% First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

%% Second-order approximation
if order_app > 1
    if linearShocksOn == 1
        mx  = nx; %[r_ba1 k_cu s_cu ]
        myx = 0; %[]
        my  = ny; 
        sigma = eta*eta';
        [D,H] = getDH_LinearShocks(n,nx,ny,nfx,nfxp,nfy,nfyp,...
            nfxx,nfxxp,nfxy,nfxyp,nfxpx,nfxpxp,nfxpy,nfxpyp,...
            nfyx,nfyxp,nfyy,nfyyp,nfypx,nfypxp,nfypy,nfypyp);
        [gxx,hxx,hss,gss] = g_h_2nd_LinearShocks(D,H,gx,hx,sigma,nx,ny,n,mx,my,myx);
    else
        % Codes provided by S. Schmitt-Grohe and M. Uribe
        [gxx,hxx] = gxx_hxx_noloop(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);
        [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
    end
else
    gxx = zeros(ny,nx,nx);
    gss = zeros(ny,1);
    hxx = zeros(nx,nx,nx);
    hss = zeros(nx,1);
end

%% Third-order approximation
if order_app > 2
    if linearShocksOn == 1
        T = getT_LinearShocks(n,nx,ny,...
            nfxxx,nfxxxp,nfxxy,nfxxyp,nfxxpxp,nfxxpy,nfxxpyp,...
            nfxyy,nfxyyp,nfxypyp,nfxpxpxp,nfxpxpy,nfxpxpyp,...
            nfxpyy,nfxpyyp,nfxpypyp,...
            nfyyy,nfyyyp,nfyypyp,nfypypyp);
        skew = zeros(nx,nx^2);
        [gxxx,hxxx,gssx,hssx,gsss,hsss] = g_h_3rd_LinearShocks(D,H,T,gx,hx,gxx,hxx,gss,hss,sigma,skew,nx,ny,n,mx,my,myx);
    else
        % The vector with the third moments
%         vectorMom3 = zeros(1,ne);
%         [gxxx,hxxx,gssx,hssx,gsss,hsss] = g_h_3rd_v3(nfx,nfxp,nfy,nfyp,...
%             nfxx,nfxxp,nfxy,nfxyp,nfxpx,nfxpxp,nfxpy,nfxpyp,...
%             nfyx,nfyxp,nfyy,nfyyp,nfypx,nfypxp,nfypy,nfypyp,...
%             nfxxx,nfxxxp,nfxxy,nfxxyp,nfxxpxp,nfxxpy,nfxxpyp,...
%             nfxyy,nfxyyp,nfxypyp,nfxpxpxp,nfxpxpy,nfxpxpyp,...
%             nfxpyy,nfxpyyp,nfxpypyp,...
%             nfyyy,nfyyyp,nfyypyp,nfypypyp,gx,gxx,gss,hx,hxx,hss,eta,vectorMom3);
        skew = zeros(nx,nx^2);
        [gxxx,hxxx,gssx,hssx,gsss,hsss] = g_h_3rd_Binning(nfx,nfxp,nfy,nfyp,...   % First-order terms
          nfxx,nfxxp,nfxy,nfxyp,nfxpx,nfxpxp,nfxpy,nfxpyp,...                                 % Second-order terms
           nfyx,nfyxp,nfyy,nfyyp,nfypx,nfypxp,nfypy,nfypyp,...
          nfxxx,nfxxxp,nfxxy,nfxxyp,nfxxpxp,nfxxpy,nfxxpyp,...                                % Third-order terms
           nfxyy,nfxyyp,nfxypyp,nfxpxpxp,nfxpxpy,nfxpxpyp,...
           nfxpyy,nfxpyyp,nfxpypyp,...
           nfyyy,nfyyyp,nfyypyp,nfypypyp,...
           gx,gxx,gss,hx,hxx,hss,eta,skew);
    end
else
    gxxx = zeros(ny,nx,nx,nx);
    gssx = zeros(ny,nx);
    gsss = zeros(ny,1);
    hxxx = zeros(nx,nx,nx,nx);    
    hssx = zeros(nx,nx);    
    hsss = zeros(nx,1);    
end

end

    

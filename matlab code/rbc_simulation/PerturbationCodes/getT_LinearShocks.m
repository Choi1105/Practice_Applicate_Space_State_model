% Copyright Andrew Binning 2013

function T = getT_LinearShocks(n,nx,ny,...
    nfxxx,nfxxxp,nfxxy,nfxxyp,nfxxpxp,nfxxpy,nfxxpyp,...
    nfxyy,nfxyyp,nfxypyp,nfxpxpxp,nfxpxpy,nfxpxpyp,...
    nfxpyy,nfxpyyp,nfxpypyp,...
    nfyyy,nfyyyp,nfyypyp,nfypypyp)

% Converts Matrices from Andreasen form to Binning form.
T = zeros(n,(2*n)^3);
nfxpypxp = permute(nfxpxpyp,[1,2,4,3]);
nfxpxxp = permute(nfxxpxp,[1,3,2,4]);
nfxpyxp = permute(nfxpxpy,[1,2,4,3]);

nfypxpxp = permute(nfxpxpyp,[1,4,2,3]);
nfypypxp = permute(nfxpypyp,[1,3,4,2]);
nfypxxp = permute(nfxxpyp,[1,4,2,3]);
nfypyxp = permute(nfxpyyp,[1,4,3,2]);

nfxxpx = permute(nfxxxp,[1,2,4,3]);
nfxypxp = permute(nfxxpyp,[1,2,4,3]);
nfxyxp = permute(nfxxpy,[1,2,4,3]);

nfyxpxp = permute(nfxpxpy,[1,4,2,3]);
nfyypxp = permute(nfxpyyp,[1,3,4,2]);
nfyxxp = permute(nfxxpy,[1,4,2,3]);
nfyyxp = permute(nfxpyy,[1,3,4,2]);

% nfxpxpyp = permute();
% nfxpypyp = permute();
nfxpxyp = permute(nfxxpyp,[1,3,2,4]);
% nfxpyyp = permute();

nfypxpyp = permute(nfxpypyp,[1,3,2,4]);
% nfypypyp = permute();
nfypxyp = permute(nfxypyp,[1,3,2,4]);
nfypyyp = permute(nfyypyp,[1,3,2,4]);

% nfxxpyp = permute();
% nfxypyp = permute();
% nfxxyp = permute();
% nfxyyp = permute();
 
nfyxpyp = permute(nfxpyyp,[1,3,2,4]);
% nfyypyp = permute();
nfyxyp = permute(nfxyyp,[1,3,2,4]);
% nfyyyp = permute();
%%%%%%%%%%%%%%
nfxpxpx = permute(nfxxpxp,[1,3,4,2]);
nfxpypx = permute(nfxxpyp,[1,3,4,2]);
nfxpxx = permute(nfxxxp,[1,4,2,3]);
nfxpyx = permute(nfxxpy,[1,3,4,2]);

nfypxpx = permute(nfxxpyp,[1,4,3,2]);
nfypypx = permute(nfxypyp,[1,4,3,2]);
nfypxx = permute(nfxxyp,[1,4,3,2]);
nfypyx = permute(nfxyyp,[1,4,3,2]);

% nfxxpx = permute();
nfxypx = permute(nfxxyp,[1,2,4,3]);
% nfxxx = permute();
nfxyx = permute(nfxxy,[1,2,4,3]);

nfyxpx = permute(nfxxpy,[1,4,3,2]);
nfyypx = permute(nfxyyp,[1,3,4,2]);
nfyxx = permute(nfxxy,[1,4,3,2]);
nfyyx = permute(nfxyy,[1,4,3,2]);
%%%%%%%%%%%%%%%%%%%%%%%%%
% nfxpxpy = permute();
nfxpypy = permute(nfxpyyp,[1,2,4,3]);
nfxpxy = permute(nfxxpy,[1,3,2,4]);
% nfxpyy = permute();

nfypxpy = permute(nfxpyyp,[1,4,2,3]);
nfypypy = permute(nfyypyp,[1,4,3,2]);
nfypxy = permute(nfxyyp,[1,4,2,3]);
nfypyy = permute(nfyyyp,[1,4,2,3]);

% nfxxpy = permute();
nfxypy = permute(nfxyyp,[1,2,4,3]);
% nfxxy = permute();
% nfxyy = permute();

nfyxpy = permute(nfxpyy,[1,3,2,4]);
nfyypy = permute(nfyyyp,[1,2,4,3]);
nfyxy = permute(nfxyy,[1,3,2,4]);
% nfyyy = permute();

for ii = 1:n
    for jj = 1:nx
        XPXP = [reshape(nfxpxpxp(ii,:,:,jj),nx,nx),reshape(nfxpypxp(ii,:,:,jj),nx,ny),reshape(nfxpxxp(ii,:,:,jj),nx,nx),reshape(nfxpyxp(ii,:,:,jj),nx,ny)]';
        XPXP = XPXP(:)';
        XPYP = [reshape(nfypxpxp(ii,:,:,jj),ny,nx),reshape(nfypypxp(ii,:,:,jj),ny,ny),reshape(nfypxxp(ii,:,:,jj),ny,nx),reshape(nfypyxp(ii,:,:,jj),ny,ny)]';
        XPYP = XPYP(:)';
        XPX = [reshape(nfxxpxp(ii,:,:,jj),nx,nx),reshape(nfxypxp(ii,:,:,jj),nx,ny),reshape(nfxxxp(ii,:,:,jj),nx,nx),reshape(nfxyxp(ii,:,:,jj),nx,ny)]';
        XPX = XPX(:)';
        XPY = [reshape(nfyxpxp(ii,:,:,jj),ny,nx),reshape(nfyypxp(ii,:,:,jj),ny,ny),reshape(nfyxxp(ii,:,:,jj),ny,nx),reshape(nfyyxp(ii,:,:,jj),ny,ny)]';
        XPY = XPY(:)';
        tmp = [XPXP,XPYP,XPX,XPY];
        T(ii,(jj-1)*(2*n)^2+(1:(2*n)^2)) = tmp;        
    end    
    for jj = 1:ny
        
        YPXP = [reshape(nfxpxpyp(ii,:,:,jj),nx,nx),reshape(nfxpypyp(ii,:,:,jj),nx,ny),reshape(nfxpxyp(ii,:,:,jj),nx,nx),reshape(nfxpyyp(ii,:,:,jj),nx,ny)]';
        YPXP = YPXP(:)';
        YPYP = [reshape(nfypxpyp(ii,:,:,jj),ny,nx),reshape(nfypypyp(ii,:,:,jj),ny,ny),reshape(nfypxyp(ii,:,:,jj),ny,nx),reshape(nfypyyp(ii,:,:,jj),ny,ny)]';
        YPYP = YPYP(:)';
        YPX = [reshape(nfxxpyp(ii,:,:,jj),nx,nx),reshape(nfxypyp(ii,:,:,jj),nx,ny),reshape(nfxxyp(ii,:,:,jj),nx,nx),reshape(nfxyyp(ii,:,:,jj),nx,ny)]';
        YPX = YPX(:)';
        YPY = [reshape(nfyxpyp(ii,:,:,jj),ny,nx),reshape(nfyypyp(ii,:,:,jj),ny,ny),reshape(nfyxyp(ii,:,:,jj),ny,nx),reshape(nfyyyp(ii,:,:,jj),ny,ny)]';
        YPY = YPY(:)';
        tmp = [YPXP,YPYP,YPX,YPY];
        T(ii,nx*(2*n)^2+(jj-1)*(2*n)^2+(1:(2*n)^2)) = tmp;
        
    end
    
    for jj = 1:nx
        XXP = [reshape(nfxpxpx(ii,:,:,jj),nx,nx),reshape(nfxpypx(ii,:,:,jj),nx,ny),reshape(nfxpxx(ii,:,:,jj),nx,nx),reshape(nfxpyx(ii,:,:,jj),nx,ny)]';
        XXP = XXP(:)';
        XYP = [reshape(nfypxpx(ii,:,:,jj),ny,nx),reshape(nfypypx(ii,:,:,jj),ny,ny),reshape(nfypxx(ii,:,:,jj),ny,nx),reshape(nfypyx(ii,:,:,jj),ny,ny)]';
        XYP = XYP(:)';
        XX = [reshape(nfxxpx(ii,:,:,jj),nx,nx),reshape(nfxypx(ii,:,:,jj),nx,ny),reshape(nfxxx(ii,:,:,jj),nx,nx),reshape(nfxyx(ii,:,:,jj),nx,ny)]';
        XX = XX(:)';
        XY = [reshape(nfyxpx(ii,:,:,jj),ny,nx),reshape(nfyypx(ii,:,:,jj),ny,ny),reshape(nfyxx(ii,:,:,jj),ny,nx),reshape(nfyyx(ii,:,:,jj),ny,ny)]';
        XY = XY(:)';
        tmp = [XXP,XYP,XX,XY];
        T(ii,n*(2*n)^2+(jj-1)*(2*n)^2+(1:(2*n)^2)) = tmp;
        
    end

   for jj = 1:ny
        YXP = [reshape(nfxpxpy(ii,:,:,jj),nx,nx),reshape(nfxpypy(ii,:,:,jj),nx,ny),reshape(nfxpxy(ii,:,:,jj),nx,nx),reshape(nfxpyy(ii,:,:,jj),nx,ny)]';
        YXP = YXP(:)';
        YYP = [reshape(nfypxpy(ii,:,:,jj),ny,nx),reshape(nfypypy(ii,:,:,jj),ny,ny),reshape(nfypxy(ii,:,:,jj),ny,nx),reshape(nfypyy(ii,:,:,jj),ny,ny)]';
        YYP = YYP(:)';
        YX = [reshape(nfxxpy(ii,:,:,jj),nx,nx),reshape(nfxypy(ii,:,:,jj),nx,ny),reshape(nfxxy(ii,:,:,jj),nx,nx),reshape(nfxyy(ii,:,:,jj),nx,ny)]';
        YX = YX(:)';
        YY = [reshape(nfyxpy(ii,:,:,jj),ny,nx),reshape(nfyypy(ii,:,:,jj),ny,ny),reshape(nfyxy(ii,:,:,jj),ny,nx),reshape(nfyyy(ii,:,:,jj),ny,ny)]';
        YY = YY(:)';
        tmp = [YXP,YYP,YX,YY];
        T(ii,(n+nx)*(2*n)^2+(jj-1)*(2*n)^2+(1:(2*n)^2)) = tmp;        
   end
end
T = sparse(T);
end
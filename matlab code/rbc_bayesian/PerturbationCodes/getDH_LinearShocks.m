% Copyright Andrew Binning 2013

function [D,H] = getDH_LinearShocks(n,nx,ny,nfx,nfxp,nfy,nfyp,... 
    nfxx,nfxxp,nfxy,nfxyp,nfxpx,nfxpxp,nfxpy,nfxpyp,...
    nfyx,nfyxp,nfyy,nfyyp,nfypx,nfypxp,nfypy,nfypyp)

% Converts Matrices from Andreasen form to Binning form.

D = [nfxp,nfyp,nfx,nfy];

H = zeros(n,(2*n)^2);
for ii = 1:n
    XP = [reshape(nfxpxp(ii,:,:,:),nx,nx),reshape(nfxpyp(ii,:,:,:),nx,ny),reshape(nfxpx(ii,:,:,:),nx,nx),reshape(nfxpy(ii,:,:,:),nx,ny)]';
    XP = XP(:)';
    YP = [reshape(nfypxp(ii,:,:,:),ny,nx),reshape(nfypyp(ii,:,:,:),ny,ny),reshape(nfypx(ii,:,:,:),ny,nx),reshape(nfypy(ii,:,:,:),ny,ny)]';
    YP = YP(:)';
    X = [reshape(nfxxp(ii,:,:,:),nx,nx),reshape(nfxyp(ii,:,:,:),nx,ny),reshape(nfxx(ii,:,:,:),nx,nx),reshape(nfxy(ii,:,:,:),nx,ny)]';
    X = X(:)';
    Y = [reshape(nfyxp(ii,:,:,:),ny,nx),reshape(nfyyp(ii,:,:,:),ny,ny),reshape(nfyx(ii,:,:,:),ny,nx),reshape(nfyy(ii,:,:,:),ny,ny)]';
    Y = Y(:)';    
    H(ii,:) = [XP,YP,X,Y];
end
H = sparse(H);
end
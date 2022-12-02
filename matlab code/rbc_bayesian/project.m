function [vcov,check] = project(vcov0,e_min,e_max)
% project -- projection of covariance matrix such the eigenvalues are
% sandwiched.
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% - **vcov0** [matrix]: initial covariance matrix
%
% - **e_min** [[]|{sqrt(eps)}]: scalar such that the minimum eigenvalue of vcov
% is greater than or equal to "e_min".
%
% - **e_max** [[]|{1/e_min}]: scalar such that maximum eigenvalue of vcov
% is less than or equal to "e_max"
%
% Outputs
% --------
%
% - **vcov** [matrix]: updated covariance matrix
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:

% if nargin<3
%     e_max=[];
%     if nargin<2
%         e_min=[];
%     end
% end

check = 0;

% if isempty(e_min)
%     e_min=sqrt(1.0e-08);
    e_min=sqrt(eps);
% end
% if isempty(e_max)
    e_max=1/e_min;
% end

vcov=.5*(vcov0+vcov0.');

% [V,D] = eig(vcov);
[V,D,W] = svd(vcov);
% d=diag(V'*vcov*V);
% D=D*spdiags(sign(d.*diag(D)),0,size(V,2),size(V,2));
D = D.*sign(diag(real(dot(V,W,1))));
% [V,D]=schur(vcov);

oldD=diag(D);

% quick exit
%------------
if any(oldD<e_min)||any(oldD>e_max)
    D=max(oldD,e_min);
    D=min(D,e_max);
    vcov = V*diag(D)*V';
    check = 1;
end


end

%Copied from Peter Acklam's website
function Z=kron3(X,Y)
%KRON   Kronecker tensor product. 
Z=reshape(permute(reshape(Y(:)*X(:).',[size(Y),size(X)]),...
[1,3,2,4]),[size(Y).*size(X)]);
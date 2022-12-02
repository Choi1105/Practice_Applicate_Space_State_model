% Copied from Van Loan
function y = KronTimesVector(B,C,x)
% y = kron(B,C)*x
[mb,nb] = size(B);
[mc,nc] = size(C);
Y = C*reshape(x,nc,nb)*B';
y = reshape(Y,mc*mb,1);
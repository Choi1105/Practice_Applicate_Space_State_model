function x=trile(A,r);

[m,n]=size(A);
[rm,rn]=size(r);
if rm~=m; error('A and r must have the same number of rows'); end;
if m~=n; error('A must be square'); end;
a(2:n)=diag(A,-1);
b=diag(A);
c(1:n-1)=diag(A,1);
if min([min(abs( [a(2:n),c(1:n-1)] )),min(abs(b))] )<1e-20;
	error('The three diagonals of A must have nonzero elements');
end;

for i=n-1:-1:1;
	cb=c(i)/b(i+1);
	b(i)=b(i)-a(i+1)*cb;
	r(i,:)=r(i,:)-r(i+1,:)*cb;
end;

A1=diag(a(2:n),-1)+diag(b);
x(1,:)=r(1,:)/b(1);
for i=2:n;
	for j=1:rn;
		x(i,j)=(r(i,j)-a(i)*x(i-1,j))/b(i);
	end;
end;

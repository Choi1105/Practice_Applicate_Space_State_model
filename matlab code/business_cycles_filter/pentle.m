function x=pentle(A,r);

[m,n]=size(A);
[rm,rn]=size(r);
if m~=n; error('A must be square'); end;
if rm~=m; error('A and r must have the same number of rows'); end;
a(3:n)=diag(A,-2);
b(2:n)=diag(A,-1);
c=diag(A);
d(1:n-1)=diag(A,1);
e(1:n-2)=diag(A,2);
if min([min(abs( [a(3:n),e(1:n-2)] )),min(abs( [b(2:n),d(1:n-1)] )),...
	min(abs(c))])<1e-20;
	error('The five diagonals of A must have nonzero elements');
end;

for i=n-2:-1:1;
	ed=e(i)/d(i+1);
	b(i)=b(i)-a(i+1)*ed;
	c(i)=c(i)-b(i+1)*ed;
	d(i)=d(i)-c(i+1)*ed;
	r(i,:)=r(i,:)-r(i+1,:)*ed;
end;

for i=3:n;
	ab=a(i)/b(i-1);
	b(i)=b(i)-c(i-1)*ab;
	c(i)=c(i)-d(i-1)*ab;
	r(i,:)=r(i,:)-r(i-1,:)*ab;
end;

A1=diag(b(2:n),-1)+diag(c)+diag(d(1:n-1),1);
x=trile(A1,r);

function [py,rflag] = priorval(theta)
% theta = cdraw;

b = [0 0.2];
c = [-0.9900    0.9900];
d = [0 0.2];
e = [-0.9900    0.9900];
f = [0 8];
g = [1 1.03];
h = [0 10];

x = zeros(length(theta),1);

x(1) = 1/(b(2)-b(1));
x(2) = 1/(c(2)-c(1));
x(3) = 1/(d(2)-d(1));
x(4) = 1/(e(2)-e(1));
x(5) = 1/(f(2)-f(1)); 
x(6) = 1/(g(2)-g(1));
x(7) = 1/(h(2)-h(1));

py = prod(x);

rflag = theta(1)<b(1) || theta(1)>b(2) ||...
    theta(2)<c(1) || theta(2)>c(2) ||...
    theta(3)<d(1) || theta(3)>d(2) ||...
    theta(4)<e(1) || theta(4)>e(2) ||...
    theta(5)<f(1) || theta(5)>f(2) ||...
    theta(6)<g(1) || theta(6)>g(2) ||...
    theta(7)<h(1) || theta(7)>h(2);

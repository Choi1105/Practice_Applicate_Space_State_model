function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = rbc.static_g1_tt(T, y, x, params);
end
g1 = zeros(12, 12);
g1(1,8)=1;
g1(1,9)=(-(params(9)+0.5*params(10)*2*(y(9)-1)));
g1(2,1)=(-1)/(y(1)*y(1))-(y(9)*y(5)+1-y(8))*params(1)*(-1)/(y(1)*y(1));
g1(2,5)=(-(y(9)*T(1)));
g1(2,8)=T(1);
g1(2,9)=(-(T(1)*y(5)));
g1(3,1)=(-1)/(y(1)*y(1))-(1+y(10))*params(1)*(-1)/(y(1)*y(1));
g1(3,10)=(-T(1));
g1(4,1)=(-(T(2)*y(11)*params(12)));
g1(4,2)=(-(y(1)*y(11)*params(12)*getPowerDeriv(y(2),params(4),1)));
g1(4,4)=1;
g1(4,11)=(-(T(2)*y(1)*params(12)));
g1(5,5)=1;
g1(5,9)=(-params(10));
g1(6,2)=(-(T(6)*T(11)));
g1(6,3)=(-(T(7)*T(4)*getPowerDeriv(y(3),params(3)-1,1)));
g1(6,5)=1;
g1(6,9)=(-(T(7)*T(5)*params(3)*y(12)*getPowerDeriv(y(9),params(3)-1,1)));
g1(6,12)=(-(T(7)*T(5)*params(3)*T(3)));
g1(7,2)=(-(y(12)*(1-params(3))*T(8)*T(9)*getPowerDeriv(y(2),(-params(3)),1)));
g1(7,3)=(-(T(10)*y(12)*(1-params(3))*T(8)*T(12)));
g1(7,4)=1;
g1(7,9)=(-(T(10)*T(9)*y(12)*(1-params(3))*T(13)));
g1(7,12)=(-(T(10)*T(9)*(1-params(3))*T(8)));
g1(8,2)=(-(T(9)*y(12)*T(8)*T(11)));
g1(8,3)=(-(T(7)*y(12)*T(8)*T(12)));
g1(8,6)=1;
g1(8,9)=(-(T(7)*T(9)*y(12)*T(13)));
g1(8,12)=(-(T(7)*T(8)*T(9)));
g1(9,1)=(-1);
g1(9,6)=1;
g1(9,7)=(-1);
g1(10,3)=1-(1-y(8));
g1(10,7)=(-1);
g1(10,8)=y(3);
g1(11,11)=1/y(11)-params(6)*1/y(11);
g1(12,12)=1/y(12)-params(5)*1/y(12);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end

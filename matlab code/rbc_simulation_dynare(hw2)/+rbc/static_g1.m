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
g1 = zeros(8, 8);
g1(1,1)=(-((-1)/(y(1)*y(1))));
g1(1,8)=1;
g1(2,2)=(-((-(T(1)*(1-params(2))*y(3)*T(6)))/((1-params(2))*y(3)*T(2)*(1-params(2))*y(3)*T(2))));
g1(2,3)=(-((-(T(1)*(1-params(2))*T(2)))/((1-params(2))*y(3)*T(2)*(1-params(2))*y(3)*T(2))));
g1(2,4)=(-(params(9)*getPowerDeriv(y(4),params(3)+params(2),1)/((1-params(2))*y(3)*T(2))));
g1(2,8)=1;
g1(3,2)=(-(T(4)*y(8)*y(6)*params(8)*params(2)*y(3)*getPowerDeriv(y(2),params(2)-1,1)));
g1(3,3)=(-(T(4)*y(8)*y(6)*params(8)*params(2)*T(5)));
g1(3,4)=(-(y(8)*y(6)*params(8)*params(2)*y(3)*T(5)*T(7)));
g1(3,6)=(-((1-params(1))*(-1)/(y(6)*y(6))+T(4)*params(2)*y(3)*T(5)*y(8)*params(8)));
g1(3,8)=1-T(4)*y(6)*params(8)*params(2)*y(3)*T(5);
g1(4,1)=(-1);
g1(4,5)=1;
g1(4,7)=(-1);
g1(5,2)=(-(T(4)*y(3)*T(6)));
g1(5,3)=(-(T(2)*T(4)));
g1(5,4)=(-(y(3)*T(2)*T(7)));
g1(5,5)=1;
g1(6,2)=(-(T(3)*(1-(1-params(1)))));
g1(6,6)=(-((y(2)-y(2)*(1-params(1)))*(-1)/(y(6)*y(6))));
g1(6,7)=1;
g1(7,3)=1/y(3)-params(4)*1/y(3);
g1(8,6)=T(3)-T(3)*params(5);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end

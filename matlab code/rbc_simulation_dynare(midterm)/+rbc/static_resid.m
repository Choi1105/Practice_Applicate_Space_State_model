function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = rbc.static_resid_tt(T, y, x, params);
end
residual = zeros(12, 1);
lhs = y(8);
rhs = params(2)+params(9)*(y(9)-1)+0.5*params(10)*(y(9)-1)^2;
residual(1) = lhs - rhs;
lhs = 1/y(1);
rhs = T(1)*(y(9)*y(5)+1-y(8));
residual(2) = lhs - rhs;
lhs = 1/y(1);
rhs = T(1)*(1+y(10));
residual(3) = lhs - rhs;
lhs = y(4);
rhs = y(1)*y(11)*params(12)*T(2);
residual(4) = lhs - rhs;
lhs = y(5);
rhs = params(9)+(y(9)-1)*params(10);
residual(5) = lhs - rhs;
lhs = y(5);
rhs = T(6)*T(7);
residual(6) = lhs - rhs;
lhs = y(4);
rhs = y(12)*(1-params(3))*T(8)*T(9)*T(10);
residual(7) = lhs - rhs;
lhs = y(6);
rhs = T(7)*T(9)*y(12)*T(8);
residual(8) = lhs - rhs;
lhs = y(6);
rhs = y(1)+y(7);
residual(9) = lhs - rhs;
lhs = y(3);
rhs = y(7)+(1-y(8))*y(3);
residual(10) = lhs - rhs;
lhs = log(y(11));
rhs = log(y(11))*params(6)+x(2);
residual(11) = lhs - rhs;
lhs = log(y(12));
rhs = log(y(12))*params(5)+x(1);
residual(12) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end

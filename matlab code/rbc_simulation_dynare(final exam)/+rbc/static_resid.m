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
residual = zeros(14, 1);
lhs = y(14);
rhs = params(8)+params(4)*(T(4)-1);
residual(1) = lhs - rhs;
lhs = y(9);
rhs = T(5)^(-params(1));
residual(2) = lhs - rhs;
lhs = y(6);
rhs = y(1)+log(y(8))+y(7);
residual(3) = lhs - rhs;
lhs = log(y(8));
rhs = y(5)-y(5)/y(14);
residual(4) = lhs - rhs;
lhs = y(6);
rhs = T(6)*y(3)*T(7);
residual(5) = lhs - rhs;
lhs = y(7);
rhs = y(2)-y(2)*(1-params(2));
residual(6) = lhs - rhs;
lhs = y(9);
rhs = y(9)*params(9)*y(14);
residual(7) = lhs - rhs;
lhs = params(6)*y(4)^(params(5)-1);
rhs = (1-params(3))*y(3)*T(8);
residual(8) = lhs - rhs;
lhs = y(9);
rhs = params(9)*y(9)*(1-params(2)+params(3)*y(3)*T(9));
residual(9) = lhs - rhs;
lhs = log(y(10));
rhs = log(y(8))/y(6);
residual(10) = lhs - rhs;
lhs = y(11);
rhs = 1;
residual(11) = lhs - rhs;
lhs = y(12);
rhs = 1;
residual(12) = lhs - rhs;
lhs = y(13);
rhs = 1;
residual(13) = lhs - rhs;
lhs = log(y(3));
rhs = log(y(3))*params(7)+x(1);
residual(14) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end

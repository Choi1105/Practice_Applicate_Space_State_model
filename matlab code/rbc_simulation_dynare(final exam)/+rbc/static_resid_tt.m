function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 9);

T(1) = ((1/params(9)-1+params(2))/params(3))^(1/(params(3)-1));
T(2) = ((1-params(3))*(T(1))^params(3)/params(6))^(1/(params(5)-1));
T(3) = (T(1)*T(2))^params(3)*(T(2))^(1-params(3));
T(4) = exp(y(5)/T(3)-(-(0.003*T(3)))/(1/params(8)-1)/T(3));
T(5) = y(1)-params(6)/params(5)*y(4)^params(5);
T(6) = y(4)^(1-params(3));
T(7) = y(2)^params(3);
T(8) = (y(2)/y(4))^params(3);
T(9) = (y(4)/y(2))^(1-params(3));

end

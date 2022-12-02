function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = rbc.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(14, 1);
lhs = y(20);
rhs = params(8)+params(4)*(T(9)-1);
residual(1) = lhs - rhs;
lhs = y(15);
rhs = T(2)^(-params(1));
residual(2) = lhs - rhs;
lhs = y(12);
rhs = y(7)+log(y(14))+y(13);
residual(3) = lhs - rhs;
lhs = log(y(14));
rhs = y(4)-y(11)/y(20);
residual(4) = lhs - rhs;
lhs = y(12);
rhs = T(3)*y(9)*T(4);
residual(5) = lhs - rhs;
lhs = y(13);
rhs = y(8)-(1-params(2))*y(2);
residual(6) = lhs - rhs;
lhs = y(15);
rhs = params(9)*y(20)*y(23);
residual(7) = lhs - rhs;
lhs = params(6)*y(10)^(params(5)-1);
rhs = (1-params(3))*y(9)*T(5);
residual(8) = lhs - rhs;
lhs = y(15);
rhs = params(9)*y(23)*(1-params(2)+params(3)*y(21)*T(6));
residual(9) = lhs - rhs;
lhs = log(y(16));
rhs = log(y(14))/y(12);
residual(10) = lhs - rhs;
lhs = y(17);
rhs = y(12)/y(5);
residual(11) = lhs - rhs;
lhs = y(18);
rhs = y(7)/y(1);
residual(12) = lhs - rhs;
lhs = y(19);
rhs = y(13)/y(6);
residual(13) = lhs - rhs;
lhs = log(y(9));
rhs = params(7)*log(y(3))+x(it_, 1);
residual(14) = lhs - rhs;

end

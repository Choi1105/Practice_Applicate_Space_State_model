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
residual = zeros(8, 1);
lhs = y(11);
rhs = 1/y(4);
residual(1) = lhs - rhs;
lhs = y(11);
rhs = T(1)/((1-params(2))*y(6)*T(5));
residual(2) = lhs - rhs;
lhs = y(11);
rhs = 1/y(14)*(1-params(1))+T(2)*y(9)*params(8)*y(15)*params(2)*y(12)*T(6);
residual(3) = lhs - rhs;
lhs = y(8);
rhs = y(4)+y(10);
residual(4) = lhs - rhs;
lhs = y(8);
rhs = T(3)*y(6)*T(5);
residual(5) = lhs - rhs;
lhs = y(10);
rhs = T(4)*(y(5)-(1-params(1))*y(1));
residual(6) = lhs - rhs;
lhs = log(y(6));
rhs = params(4)*log(y(2))+x(it_, 1);
residual(7) = lhs - rhs;
lhs = log(y(9));
rhs = params(5)*log(y(3))+x(it_, 2);
residual(8) = lhs - rhs;

end

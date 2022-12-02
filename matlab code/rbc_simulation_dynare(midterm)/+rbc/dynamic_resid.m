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
residual = zeros(12, 1);
lhs = y(13);
rhs = params(2)+params(9)*(y(14)-1)+0.5*params(10)*(y(14)-1)^2;
residual(1) = lhs - rhs;
lhs = 1/y(6);
rhs = T(1)*(y(19)*y(21)+1-y(20));
residual(2) = lhs - rhs;
lhs = 1/y(6);
rhs = T(1)*(1+y(15));
residual(3) = lhs - rhs;
lhs = y(9);
rhs = y(6)*y(16)*params(12)*T(2);
residual(4) = lhs - rhs;
lhs = y(10);
rhs = params(9)+(y(14)-1)*params(10);
residual(5) = lhs - rhs;
lhs = y(10);
rhs = T(6)*T(7);
residual(6) = lhs - rhs;
lhs = y(9);
rhs = y(17)*(1-params(3))*T(8)*T(9)*T(10);
residual(7) = lhs - rhs;
lhs = y(11);
rhs = T(7)*T(9)*y(17)*T(8);
residual(8) = lhs - rhs;
lhs = y(11);
rhs = y(6)+y(12);
residual(9) = lhs - rhs;
lhs = y(8);
rhs = y(2)+(1-y(3))*y(1);
residual(10) = lhs - rhs;
lhs = log(y(16));
rhs = params(6)*log(y(4))+x(it_, 2);
residual(11) = lhs - rhs;
lhs = log(y(17));
rhs = params(5)*log(y(5))+x(it_, 1);
residual(12) = lhs - rhs;

end

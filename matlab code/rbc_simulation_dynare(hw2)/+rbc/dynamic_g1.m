function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = rbc.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(8, 17);
g1(1,4)=(-((-1)/(y(4)*y(4))));
g1(1,11)=1;
g1(2,1)=(-((-(T(1)*(1-params(2))*y(6)*T(7)))/((1-params(2))*y(6)*T(5)*(1-params(2))*y(6)*T(5))));
g1(2,6)=(-((-(T(1)*(1-params(2))*T(5)))/((1-params(2))*y(6)*T(5)*(1-params(2))*y(6)*T(5))));
g1(2,7)=(-(params(9)*getPowerDeriv(y(7),params(3)+params(2),1)/((1-params(2))*y(6)*T(5))));
g1(2,11)=1;
g1(3,5)=(-(T(2)*y(9)*params(8)*y(15)*params(2)*y(12)*getPowerDeriv(y(5),params(2)-1,1)));
g1(3,12)=(-(T(2)*y(9)*params(8)*y(15)*params(2)*T(6)));
g1(3,13)=(-(y(9)*params(8)*y(15)*params(2)*y(12)*T(6)*getPowerDeriv(y(13),1-params(2),1)));
g1(3,9)=(-(T(2)*params(2)*y(12)*T(6)*params(8)*y(15)));
g1(3,14)=(-((1-params(1))*(-1)/(y(14)*y(14))));
g1(3,11)=1;
g1(3,15)=(-(T(2)*y(9)*params(8)*params(2)*y(12)*T(6)));
g1(4,4)=(-1);
g1(4,8)=1;
g1(4,10)=(-1);
g1(5,1)=(-(T(3)*y(6)*T(7)));
g1(5,6)=(-(T(3)*T(5)));
g1(5,7)=(-(y(6)*T(5)*getPowerDeriv(y(7),1-params(2),1)));
g1(5,8)=1;
g1(6,1)=(-(T(4)*(-(1-params(1)))));
g1(6,5)=(-T(4));
g1(6,9)=(-((y(5)-(1-params(1))*y(1))*(-1)/(y(9)*y(9))));
g1(6,10)=1;
g1(7,2)=(-(params(4)*1/y(2)));
g1(7,6)=1/y(6);
g1(7,16)=(-1);
g1(8,3)=(-(params(5)*1/y(3)));
g1(8,9)=T(4);
g1(8,17)=(-1);

end

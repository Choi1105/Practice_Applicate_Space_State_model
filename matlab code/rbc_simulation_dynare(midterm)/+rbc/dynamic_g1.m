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
g1 = zeros(12, 23);
g1(1,13)=1;
g1(1,14)=(-(params(9)+0.5*params(10)*2*(y(14)-1)));
g1(2,6)=(-1)/(y(6)*y(6));
g1(2,18)=(-((y(19)*y(21)+1-y(20))*params(1)*(-1)/(y(18)*y(18))));
g1(2,19)=(-(T(1)*y(21)));
g1(2,20)=T(1);
g1(2,21)=(-(T(1)*y(19)));
g1(3,6)=(-1)/(y(6)*y(6));
g1(3,18)=(-((1+y(15))*params(1)*(-1)/(y(18)*y(18))));
g1(3,15)=(-T(1));
g1(4,6)=(-(T(2)*y(16)*params(12)));
g1(4,7)=(-(y(6)*y(16)*params(12)*getPowerDeriv(y(7),params(4),1)));
g1(4,9)=1;
g1(4,16)=(-(T(2)*y(6)*params(12)));
g1(5,10)=1;
g1(5,14)=(-params(10));
g1(6,7)=(-(T(6)*T(11)));
g1(6,8)=(-(T(7)*T(4)*getPowerDeriv(y(8),params(3)-1,1)));
g1(6,10)=1;
g1(6,14)=(-(T(7)*T(5)*params(3)*y(17)*getPowerDeriv(y(14),params(3)-1,1)));
g1(6,17)=(-(T(7)*T(5)*params(3)*T(3)));
g1(7,7)=(-(y(17)*(1-params(3))*T(8)*T(9)*getPowerDeriv(y(7),(-params(3)),1)));
g1(7,8)=(-(T(10)*y(17)*(1-params(3))*T(8)*T(12)));
g1(7,9)=1;
g1(7,14)=(-(T(10)*T(9)*y(17)*(1-params(3))*T(13)));
g1(7,17)=(-(T(10)*T(9)*(1-params(3))*T(8)));
g1(8,7)=(-(T(9)*y(17)*T(8)*T(11)));
g1(8,8)=(-(T(7)*y(17)*T(8)*T(12)));
g1(8,11)=1;
g1(8,14)=(-(T(7)*T(9)*y(17)*T(13)));
g1(8,17)=(-(T(7)*T(8)*T(9)));
g1(9,6)=(-1);
g1(9,11)=1;
g1(9,12)=(-1);
g1(10,1)=(-(1-y(3)));
g1(10,8)=1;
g1(10,2)=(-1);
g1(10,3)=y(1);
g1(11,4)=(-(params(6)*1/y(4)));
g1(11,16)=1/y(16);
g1(11,23)=(-1);
g1(12,5)=(-(params(5)*1/y(5)));
g1(12,17)=1/y(17);
g1(12,22)=(-1);

end

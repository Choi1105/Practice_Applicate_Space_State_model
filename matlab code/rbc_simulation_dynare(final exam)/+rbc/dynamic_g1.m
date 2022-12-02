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
g1 = zeros(14, 24);
g1(1,4)=(-(params(4)*T(9)*1/T(8)));
g1(1,20)=1;
g1(2,7)=(-T(10));
g1(2,10)=(-(T(10)*(-(params(6)/params(5)*getPowerDeriv(y(10),params(5),1)))));
g1(2,15)=1;
g1(3,7)=(-1);
g1(3,12)=1;
g1(3,13)=(-1);
g1(3,14)=(-(1/y(14)));
g1(4,4)=(-1);
g1(4,11)=1/y(20);
g1(4,14)=1/y(14);
g1(4,20)=(-y(11))/(y(20)*y(20));
g1(5,2)=(-(T(3)*y(9)*getPowerDeriv(y(2),params(3),1)));
g1(5,9)=(-(T(3)*T(4)));
g1(5,10)=(-(y(9)*T(4)*getPowerDeriv(y(10),1-params(3),1)));
g1(5,12)=1;
g1(6,2)=1-params(2);
g1(6,8)=(-1);
g1(6,13)=1;
g1(7,15)=1;
g1(7,23)=(-(params(9)*y(20)));
g1(7,20)=(-(params(9)*y(23)));
g1(8,2)=(-((1-params(3))*y(9)*1/y(10)*T(11)));
g1(8,9)=(-((1-params(3))*T(5)));
g1(8,10)=params(6)*getPowerDeriv(y(10),params(5)-1,1)-(1-params(3))*y(9)*T(11)*(-y(2))/(y(10)*y(10));
g1(9,8)=(-(params(9)*y(23)*params(3)*y(21)*(-y(22))/(y(8)*y(8))*T(12)));
g1(9,21)=(-(params(9)*y(23)*params(3)*T(6)));
g1(9,22)=(-(params(9)*y(23)*params(3)*y(21)*T(12)*1/y(8)));
g1(9,15)=1;
g1(9,23)=(-(params(9)*(1-params(2)+params(3)*y(21)*T(6))));
g1(10,12)=(-((-log(y(14)))/(y(12)*y(12))));
g1(10,14)=(-(1/y(14)/y(12)));
g1(10,16)=1/y(16);
g1(11,5)=(-((-y(12))/(y(5)*y(5))));
g1(11,12)=(-(1/y(5)));
g1(11,17)=1;
g1(12,1)=(-((-y(7))/(y(1)*y(1))));
g1(12,7)=(-(1/y(1)));
g1(12,18)=1;
g1(13,6)=(-((-y(13))/(y(6)*y(6))));
g1(13,13)=(-(1/y(6)));
g1(13,19)=1;
g1(14,3)=(-(params(7)*1/y(3)));
g1(14,9)=1/y(9);
g1(14,24)=(-1);

end

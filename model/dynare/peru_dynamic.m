function [residual, g1, g2, g3] = peru_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(33, 1);
T26 = (1-params(2))/params(2)*(params(3)-params(3)*params(4))/(params(5)/params(6))*params(7)/params(8);
T129 = params(14)*(params(8)/params(15))^(params(16)-1);
T141 = params(14)*(params(8)*params(17)/params(15))^(params(16)-1);
T181 = params(41)/(params(15)*params(38));
T202 = 1/(params(15)*params(38));
lhs =y(22)-y(23);
rhs =T26*((y(24)-params(4)*y(2))/(1-params(4))-y(25))+params(1)*(1-params(6))*((y(24)-params(4)*y(2)-(y(57)-params(4)*y(24)))/(1-params(4))+y(55)-y(56));
residual(1)= lhs-rhs;
lhs =params(4)*(y(24)-y(2));
rhs =y(57)-y(24)-(1-params(4))*(y(26)-(y(56)-y(23)));
residual(2)= lhs-rhs;
lhs =params(9)*y(27);
rhs =y(28)-y(23)-(y(24)-params(4)*y(2))/(1-params(4));
residual(3)= lhs-rhs;
lhs =y(29);
rhs =y(27)+y(30);
residual(4)= lhs-rhs;
lhs =y(43)-y(3);
rhs =params(10)*(y(31)-y(3));
residual(5)= lhs-rhs;
lhs =y(32);
rhs =(1-params(11))/params(11)*(y(33)-y(34));
residual(6)= lhs-rhs;
lhs =y(33);
rhs =(1-params(1)*params(11))*(y(29)+y(35))+params(1)*params(11)*(params(12)*y(58)+y(59));
residual(7)= lhs-rhs;
lhs =y(34);
rhs =y(29)*(1-params(1)*params(11))+params(1)*params(11)*(y(58)*(params(12)-1)+y(60));
residual(8)= lhs-rhs;
lhs =y(35);
rhs =y(28)-y(36)-y(30);
residual(9)= lhs-rhs;
lhs =y(22)-y(23);
rhs =(y(31)-y(3))*params(13);
residual(10)= lhs-rhs;
lhs =y(41);
rhs =y(36)-y(37)-y(23);
residual(11)= lhs-rhs;
lhs =y(36)-y(23);
rhs =y(37)*T129;
residual(12)= lhs-rhs;
lhs =y(44);
rhs =y(32)+y(7)-y(1)-y(37)*T141;
residual(13)= lhs-rhs;
lhs =y(37);
rhs =y(32)+y(8)-y(45)-y(38)+y(9);
residual(14)= lhs-rhs;
lhs =y(29);
rhs =(y(36)-y(23))*(-params(16))+((1-params(14))*(params(3)*y(24)+y(31)*params(18))+params(14)*params(19)*params(20)*(y(41)*params(16)+y(46)))/(params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)));
residual(15)= lhs-rhs;
lhs =y(29)+y(36)-y(23);
rhs =T181*(y(23)+params(21)*(y(10)+y(4)-y(23))-y(39)-params(21)*(y(26)+y(39)-y(23))+params(22)*(y(39)+y(40)-y(23)))+T202*(y(24)*params(3)*params(8)+params(5)*params(7)*(y(22)-y(23)+y(43)));
residual(16)= lhs-rhs;
lhs =y(40);
rhs =y(47)+y(61)-y(38);
residual(17)= lhs-rhs;
lhs =y(26);
rhs =y(40)+y(48)+params(23)*(y(39)-y(22)-y(25))+y(42);
residual(18)= lhs-rhs;
lhs =y(42);
rhs =params(24)*(y(10)-y(16));
residual(19)= lhs-rhs;
lhs =y(40);
rhs =params(25)*y(11)+(1-params(25))*(params(26)*(y(1)-y(17))+params(27)*(y(5)-y(18))+params(28)*(y(10)-y(16)));
residual(20)= lhs-rhs;
lhs =y(44);
rhs =y(23)-y(1);
residual(21)= lhs-rhs;
lhs =y(43);
rhs =1/params(6)*(y(25)-(1-params(6))*y(3));
residual(22)= lhs-rhs;
lhs =y(30);
rhs =x(it_, 1)+params(29)*y(6)+params(30)*y(19);
residual(23)= lhs-rhs;
lhs =y(48);
rhs =params(31)*y(15)+x(it_, 2);
residual(24)= lhs-rhs;
lhs =y(45);
rhs =params(32)*y(12)+x(it_, 3);
residual(25)= lhs-rhs;
lhs =y(46);
rhs =x(it_, 4)+params(33)*y(13)+params(35)*y(20);
residual(26)= lhs-rhs;
lhs =y(47);
rhs =x(it_, 5)+params(34)*y(14)+params(36)*y(21);
residual(27)= lhs-rhs;
lhs =y(49);
rhs =y(10);
residual(28)= lhs-rhs;
lhs =y(50);
rhs =y(1);
residual(29)= lhs-rhs;
lhs =y(51);
rhs =y(5);
residual(30)= lhs-rhs;
lhs =y(52);
rhs =y(6);
residual(31)= lhs-rhs;
lhs =y(53);
rhs =y(13);
residual(32)= lhs-rhs;
lhs =y(54);
rhs =y(14);
residual(33)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(33, 66);

  %
  % Jacobian matrix
  %

  g1(1,22)=1;
  g1(1,55)=(-(params(1)*(1-params(6))));
  g1(1,23)=(-1);
  g1(1,56)=params(1)*(1-params(6));
  g1(1,2)=(-(T26*(-params(4))/(1-params(4))+params(1)*(1-params(6))*(-params(4))/(1-params(4))));
  g1(1,24)=(-(T26*1/(1-params(4))+params(1)*(1-params(6))*(1-(-params(4)))/(1-params(4))));
  g1(1,57)=(-(params(1)*(1-params(6))*(-1)/(1-params(4))));
  g1(1,25)=T26;
  g1(2,23)=1-params(4);
  g1(2,56)=(-(1-params(4)));
  g1(2,2)=(-params(4));
  g1(2,24)=params(4)-(-1);
  g1(2,57)=(-1);
  g1(2,26)=1-params(4);
  g1(3,23)=1;
  g1(3,2)=(-params(4))/(1-params(4));
  g1(3,24)=1/(1-params(4));
  g1(3,27)=params(9);
  g1(3,28)=(-1);
  g1(4,27)=(-1);
  g1(4,29)=1;
  g1(4,30)=(-1);
  g1(5,3)=(-1)-(-params(10));
  g1(5,31)=(-params(10));
  g1(5,43)=1;
  g1(6,32)=1;
  g1(6,33)=(-((1-params(11))/params(11)));
  g1(6,34)=(1-params(11))/params(11);
  g1(7,29)=(-(1-params(1)*params(11)));
  g1(7,58)=(-(params(1)*params(11)*params(12)));
  g1(7,33)=1;
  g1(7,59)=(-(params(1)*params(11)));
  g1(7,35)=(-(1-params(1)*params(11)));
  g1(8,29)=(-(1-params(1)*params(11)));
  g1(8,58)=(-(params(1)*params(11)*(params(12)-1)));
  g1(8,34)=1;
  g1(8,60)=(-(params(1)*params(11)));
  g1(9,28)=(-1);
  g1(9,30)=1;
  g1(9,35)=1;
  g1(9,36)=1;
  g1(10,22)=1;
  g1(10,23)=(-1);
  g1(10,3)=params(13);
  g1(10,31)=(-params(13));
  g1(11,23)=1;
  g1(11,36)=(-1);
  g1(11,37)=1;
  g1(11,41)=1;
  g1(12,23)=(-1);
  g1(12,36)=1;
  g1(12,37)=(-T129);
  g1(13,1)=1;
  g1(13,32)=(-1);
  g1(13,7)=(-1);
  g1(13,37)=T141;
  g1(13,44)=1;
  g1(14,32)=(-1);
  g1(14,8)=(-1);
  g1(14,37)=1;
  g1(14,9)=(-1);
  g1(14,38)=1;
  g1(14,45)=1;
  g1(15,23)=(-params(16));
  g1(15,24)=(-(params(3)*(1-params(14))/(params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,29)=1;
  g1(15,31)=(-((1-params(14))*params(18)/(params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,36)=params(16);
  g1(15,41)=(-(params(16)*params(14)*params(19)*params(20)/(params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,46)=(-(params(14)*params(19)*params(20)/(params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(16,22)=(-(T202*params(5)*params(7)));
  g1(16,23)=(-1)-(T181*(1-params(21)-(-params(21))-params(22))+T202*(-(params(5)*params(7))));
  g1(16,24)=(-(T202*params(3)*params(8)));
  g1(16,4)=(-(T181*params(21)));
  g1(16,26)=(-(T181*(-params(21))));
  g1(16,29)=1;
  g1(16,36)=1;
  g1(16,10)=(-(T181*params(21)));
  g1(16,39)=(-(T181*(params(22)+(-1)-params(21))));
  g1(16,40)=(-(T181*params(22)));
  g1(16,43)=(-(T202*params(5)*params(7)));
  g1(17,38)=1;
  g1(17,61)=(-1);
  g1(17,40)=1;
  g1(17,47)=(-1);
  g1(18,22)=params(23);
  g1(18,25)=params(23);
  g1(18,26)=1;
  g1(18,39)=(-params(23));
  g1(18,40)=(-1);
  g1(18,42)=(-1);
  g1(18,48)=(-1);
  g1(19,10)=(-params(24));
  g1(19,42)=1;
  g1(19,16)=params(24);
  g1(20,1)=(-((1-params(25))*params(26)));
  g1(20,5)=(-((1-params(25))*params(27)));
  g1(20,10)=(-((1-params(25))*params(28)));
  g1(20,11)=(-params(25));
  g1(20,40)=1;
  g1(20,16)=(-((1-params(25))*(-params(28))));
  g1(20,17)=(-((1-params(25))*(-params(26))));
  g1(20,18)=(-((1-params(25))*(-params(27))));
  g1(21,1)=1;
  g1(21,23)=(-1);
  g1(21,44)=1;
  g1(22,3)=(-(1/params(6)*(-(1-params(6)))));
  g1(22,25)=(-(1/params(6)));
  g1(22,43)=1;
  g1(23,6)=(-params(29));
  g1(23,30)=1;
  g1(23,62)=(-1);
  g1(23,19)=(-params(30));
  g1(24,15)=(-params(31));
  g1(24,48)=1;
  g1(24,63)=(-1);
  g1(25,12)=(-params(32));
  g1(25,45)=1;
  g1(25,64)=(-1);
  g1(26,13)=(-params(33));
  g1(26,46)=1;
  g1(26,65)=(-1);
  g1(26,20)=(-params(35));
  g1(27,14)=(-params(34));
  g1(27,47)=1;
  g1(27,66)=(-1);
  g1(27,21)=(-params(36));
  g1(28,10)=(-1);
  g1(28,49)=1;
  g1(29,1)=(-1);
  g1(29,50)=1;
  g1(30,5)=(-1);
  g1(30,51)=1;
  g1(31,6)=(-1);
  g1(31,52)=1;
  g1(32,13)=(-1);
  g1(32,53)=1;
  g1(33,14)=(-1);
  g1(33,54)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],33,4356);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],33,287496);
end
end
end
end

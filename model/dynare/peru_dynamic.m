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

residual = zeros(30, 1);
T26 = (1-params(2))/params(2)*(params(3)-params(3)*params(4))/(params(5)/params(6))*params(7)/params(8);
T129 = params(14)*(params(8)/params(15))^(params(16)-1);
T141 = params(14)*(params(8)*params(17)/params(15))^(params(16)-1);
T182 = params(41)/(params(15)*params(38));
T203 = 1/(params(15)*params(38));
lhs =y(19)-y(20);
rhs =T26*((y(21)-params(4)*y(2))/(1-params(4))-y(22))+params(1)*(1-params(6))*((y(21)-params(4)*y(2)-(y(51)-params(4)*y(21)))/(1-params(4))+y(49)-y(50));
residual(1)= lhs-rhs;
lhs =params(4)*(y(21)-y(2));
rhs =y(51)-y(21)-(1-params(4))*(y(23)-(y(50)-y(20)));
residual(2)= lhs-rhs;
lhs =params(9)*y(24);
rhs =y(25)-y(20)-(y(21)-params(4)*y(2))/(1-params(4));
residual(3)= lhs-rhs;
lhs =y(26);
rhs =y(24)+y(27);
residual(4)= lhs-rhs;
lhs =y(40)-y(3);
rhs =params(10)*(y(28)-y(3));
residual(5)= lhs-rhs;
lhs =y(29);
rhs =(1-params(11))/params(11)*(y(30)-y(31));
residual(6)= lhs-rhs;
lhs =y(30);
rhs =(1-params(1)*params(11))*(y(26)+y(32))+params(1)*params(11)*(params(12)*y(52)+y(53));
residual(7)= lhs-rhs;
lhs =y(31);
rhs =y(26)*(1-params(1)*params(11))+params(1)*params(11)*(y(52)*(params(12)-1)+y(54));
residual(8)= lhs-rhs;
lhs =y(32);
rhs =y(25)-y(33)-y(27);
residual(9)= lhs-rhs;
lhs =y(19)-y(20);
rhs =(y(28)-y(3))*params(13);
residual(10)= lhs-rhs;
lhs =y(38);
rhs =y(33)-y(34)-y(20);
residual(11)= lhs-rhs;
lhs =y(33)-y(20);
rhs =y(34)*T129;
residual(12)= lhs-rhs;
lhs =y(41);
rhs =y(29)+y(7)-y(1)-y(34)*T141;
residual(13)= lhs-rhs;
lhs =y(34);
rhs =y(29)+y(8)-y(42)-y(35)+y(9);
residual(14)= lhs-rhs;
lhs =y(26);
rhs =(y(33)-y(20))*(-params(16))+((1-params(14))*(params(3)*y(21)+y(28)*params(18))+params(16)*params(14)*params(19)*params(20)*(y(38)*params(16)+y(43)))/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)));
residual(15)= lhs-rhs;
lhs =y(26)+y(33)-y(20);
rhs =T182*(y(20)+params(21)*(y(10)+y(4)-y(20))-y(36)-params(21)*(y(23)+y(36)-y(20))+params(22)*(y(36)+y(37)-y(20)))+T203*(y(21)*params(3)*params(8)+params(5)*params(7)*(y(19)-y(20)+y(40)));
residual(16)= lhs-rhs;
lhs =y(37);
rhs =y(44)+y(55)-y(35);
residual(17)= lhs-rhs;
lhs =y(23);
rhs =y(37)+y(45)+params(23)*(y(36)-y(19)-y(22))+y(39);
residual(18)= lhs-rhs;
lhs =y(39);
rhs =params(24)*(y(10)-y(16));
residual(19)= lhs-rhs;
lhs =y(37);
rhs =params(25)*y(11)+(1-params(25))*(params(26)*(y(1)-y(17))+params(27)*(y(5)-y(18))+params(28)*(y(10)-y(16)));
residual(20)= lhs-rhs;
lhs =y(41);
rhs =y(20)-y(1);
residual(21)= lhs-rhs;
lhs =y(40);
rhs =1/params(6)*(y(22)-(1-params(6))*y(3));
residual(22)= lhs-rhs;
lhs =y(27);
rhs =params(29)*y(6)+x(it_, 1);
residual(23)= lhs-rhs;
lhs =y(45);
rhs =params(31)*y(15)+x(it_, 2);
residual(24)= lhs-rhs;
lhs =y(42);
rhs =params(32)*y(12)+x(it_, 3);
residual(25)= lhs-rhs;
lhs =y(43);
rhs =params(33)*y(13)+x(it_, 4);
residual(26)= lhs-rhs;
lhs =y(44);
rhs =params(34)*y(14)+x(it_, 5);
residual(27)= lhs-rhs;
lhs =y(46);
rhs =y(10);
residual(28)= lhs-rhs;
lhs =y(47);
rhs =y(1);
residual(29)= lhs-rhs;
lhs =y(48);
rhs =y(5);
residual(30)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(30, 60);

  %
  % Jacobian matrix
  %

  g1(1,19)=1;
  g1(1,49)=(-(params(1)*(1-params(6))));
  g1(1,20)=(-1);
  g1(1,50)=params(1)*(1-params(6));
  g1(1,2)=(-(T26*(-params(4))/(1-params(4))+params(1)*(1-params(6))*(-params(4))/(1-params(4))));
  g1(1,21)=(-(T26*1/(1-params(4))+params(1)*(1-params(6))*(1-(-params(4)))/(1-params(4))));
  g1(1,51)=(-(params(1)*(1-params(6))*(-1)/(1-params(4))));
  g1(1,22)=T26;
  g1(2,20)=1-params(4);
  g1(2,50)=(-(1-params(4)));
  g1(2,2)=(-params(4));
  g1(2,21)=params(4)-(-1);
  g1(2,51)=(-1);
  g1(2,23)=1-params(4);
  g1(3,20)=1;
  g1(3,2)=(-params(4))/(1-params(4));
  g1(3,21)=1/(1-params(4));
  g1(3,24)=params(9);
  g1(3,25)=(-1);
  g1(4,24)=(-1);
  g1(4,26)=1;
  g1(4,27)=(-1);
  g1(5,3)=(-1)-(-params(10));
  g1(5,28)=(-params(10));
  g1(5,40)=1;
  g1(6,29)=1;
  g1(6,30)=(-((1-params(11))/params(11)));
  g1(6,31)=(1-params(11))/params(11);
  g1(7,26)=(-(1-params(1)*params(11)));
  g1(7,52)=(-(params(1)*params(11)*params(12)));
  g1(7,30)=1;
  g1(7,53)=(-(params(1)*params(11)));
  g1(7,32)=(-(1-params(1)*params(11)));
  g1(8,26)=(-(1-params(1)*params(11)));
  g1(8,52)=(-(params(1)*params(11)*(params(12)-1)));
  g1(8,31)=1;
  g1(8,54)=(-(params(1)*params(11)));
  g1(9,25)=(-1);
  g1(9,27)=1;
  g1(9,32)=1;
  g1(9,33)=1;
  g1(10,19)=1;
  g1(10,20)=(-1);
  g1(10,3)=params(13);
  g1(10,28)=(-params(13));
  g1(11,20)=1;
  g1(11,33)=(-1);
  g1(11,34)=1;
  g1(11,38)=1;
  g1(12,20)=(-1);
  g1(12,33)=1;
  g1(12,34)=(-T129);
  g1(13,1)=1;
  g1(13,29)=(-1);
  g1(13,7)=(-1);
  g1(13,34)=T141;
  g1(13,41)=1;
  g1(14,29)=(-1);
  g1(14,8)=(-1);
  g1(14,34)=1;
  g1(14,9)=(-1);
  g1(14,35)=1;
  g1(14,42)=1;
  g1(15,20)=(-params(16));
  g1(15,21)=(-(params(3)*(1-params(14))/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,26)=1;
  g1(15,28)=(-((1-params(14))*params(18)/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,33)=params(16);
  g1(15,38)=(-(params(16)*params(16)*params(14)*params(19)*params(20)/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,43)=(-(params(16)*params(14)*params(19)*params(20)/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(16,19)=(-(T203*params(5)*params(7)));
  g1(16,20)=(-1)-(T182*(1-params(21)-(-params(21))-params(22))+T203*(-(params(5)*params(7))));
  g1(16,21)=(-(T203*params(3)*params(8)));
  g1(16,4)=(-(T182*params(21)));
  g1(16,23)=(-(T182*(-params(21))));
  g1(16,26)=1;
  g1(16,33)=1;
  g1(16,10)=(-(T182*params(21)));
  g1(16,36)=(-(T182*(params(22)+(-1)-params(21))));
  g1(16,37)=(-(T182*params(22)));
  g1(16,40)=(-(T203*params(5)*params(7)));
  g1(17,35)=1;
  g1(17,55)=(-1);
  g1(17,37)=1;
  g1(17,44)=(-1);
  g1(18,19)=params(23);
  g1(18,22)=params(23);
  g1(18,23)=1;
  g1(18,36)=(-params(23));
  g1(18,37)=(-1);
  g1(18,39)=(-1);
  g1(18,45)=(-1);
  g1(19,10)=(-params(24));
  g1(19,39)=1;
  g1(19,16)=params(24);
  g1(20,1)=(-((1-params(25))*params(26)));
  g1(20,5)=(-((1-params(25))*params(27)));
  g1(20,10)=(-((1-params(25))*params(28)));
  g1(20,11)=(-params(25));
  g1(20,37)=1;
  g1(20,16)=(-((1-params(25))*(-params(28))));
  g1(20,17)=(-((1-params(25))*(-params(26))));
  g1(20,18)=(-((1-params(25))*(-params(27))));
  g1(21,1)=1;
  g1(21,20)=(-1);
  g1(21,41)=1;
  g1(22,3)=(-(1/params(6)*(-(1-params(6)))));
  g1(22,22)=(-(1/params(6)));
  g1(22,40)=1;
  g1(23,6)=(-params(29));
  g1(23,27)=1;
  g1(23,56)=(-1);
  g1(24,15)=(-params(31));
  g1(24,45)=1;
  g1(24,57)=(-1);
  g1(25,12)=(-params(32));
  g1(25,42)=1;
  g1(25,58)=(-1);
  g1(26,13)=(-params(33));
  g1(26,43)=1;
  g1(26,59)=(-1);
  g1(27,14)=(-params(34));
  g1(27,44)=1;
  g1(27,60)=(-1);
  g1(28,10)=(-1);
  g1(28,46)=1;
  g1(29,1)=(-1);
  g1(29,47)=1;
  g1(30,5)=(-1);
  g1(30,48)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],30,3600);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],30,216000);
end
end
end
end

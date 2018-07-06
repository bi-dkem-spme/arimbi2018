function [residual, g1, g2, g3] = peru_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 30, 1);

%
% Model equations
%

T26 = (1-params(2))/params(2)*(params(3)-params(3)*params(4))/(params(5)/params(6))*params(7)/params(8);
T109 = params(14)*(params(8)/params(15))^(params(16)-1);
T119 = params(14)*(params(8)*params(17)/params(15))^(params(16)-1);
T158 = params(41)/(params(15)*params(38));
T174 = 1/(params(15)*params(38));
lhs =y(1)-y(2);
rhs =T26*((y(3)-params(4)*y(3))/(1-params(4))-y(4))+(y(1)-y(2))*params(1)*(1-params(6));
residual(1)= lhs-rhs;
lhs =0;
rhs =(-((1-params(4))*y(5)));
residual(2)= lhs-rhs;
lhs =params(9)*y(6);
rhs =y(7)-y(2)-(y(3)-params(4)*y(3))/(1-params(4));
residual(3)= lhs-rhs;
lhs =y(8);
rhs =y(6)+y(9);
residual(4)= lhs-rhs;
lhs =y(22)-y(4);
rhs =params(10)*(y(10)-y(4));
residual(5)= lhs-rhs;
lhs =y(11);
rhs =(1-params(11))/params(11)*(y(12)-y(13));
residual(6)= lhs-rhs;
lhs =y(12);
rhs =(1-params(1)*params(11))*(y(8)+y(14))+params(1)*params(11)*(y(12)+y(11)*params(12));
residual(7)= lhs-rhs;
lhs =y(13);
rhs =y(8)*(1-params(1)*params(11))+params(1)*params(11)*(y(13)+y(11)*(params(12)-1));
residual(8)= lhs-rhs;
lhs =y(14);
rhs =y(7)-y(15)-y(9);
residual(9)= lhs-rhs;
lhs =y(1)-y(2);
rhs =(y(10)-y(4))*params(13);
residual(10)= lhs-rhs;
lhs =y(20);
rhs =y(15)-y(16)-y(2);
residual(11)= lhs-rhs;
lhs =y(15)-y(2);
rhs =y(16)*T109;
residual(12)= lhs-rhs;
lhs =y(23);
rhs =y(11)+y(15)-y(2)-y(16)*T119;
residual(13)= lhs-rhs;
lhs =y(16);
rhs =y(17)+y(11)+y(16)-y(24)-y(17);
residual(14)= lhs-rhs;
lhs =y(8);
rhs =(y(15)-y(2))*(-params(16))+((1-params(14))*(params(3)*y(3)+y(10)*params(18))+params(16)*params(14)*params(19)*params(20)*(y(20)*params(16)+y(25)))/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)));
residual(15)= lhs-rhs;
lhs =y(8)+y(15)-y(2);
rhs =T158*(y(2)+params(21)*(y(5)+y(18)-y(2))-y(18)-params(21)*(y(5)+y(18)-y(2))+params(22)*(y(18)+y(19)-y(2)))+T174*(y(3)*params(3)*params(8)+params(5)*params(7)*(y(1)-y(2)+y(22)));
residual(16)= lhs-rhs;
lhs =y(19);
rhs =y(17)+y(26)-y(17);
residual(17)= lhs-rhs;
lhs =y(5);
rhs =y(19)+y(27)+params(23)*(y(18)-y(1)-y(4))+y(21);
residual(18)= lhs-rhs;
residual(19) = y(21);
lhs =y(19);
rhs =y(19)*params(25);
residual(20)= lhs-rhs;
residual(21) = y(23);
lhs =y(22);
rhs =1/params(6)*(y(4)-y(4)*(1-params(6)));
residual(22)= lhs-rhs;
lhs =y(9);
rhs =y(9)*params(29)+x(1);
residual(23)= lhs-rhs;
lhs =y(27);
rhs =y(27)*params(31)+x(2);
residual(24)= lhs-rhs;
lhs =y(24);
rhs =y(24)*params(32)+x(3);
residual(25)= lhs-rhs;
lhs =y(25);
rhs =y(25)*params(33)+x(4);
residual(26)= lhs-rhs;
lhs =y(26);
rhs =y(26)*params(34)+x(5);
residual(27)= lhs-rhs;
lhs =y(28);
rhs =y(18);
residual(28)= lhs-rhs;
lhs =y(29);
rhs =y(2);
residual(29)= lhs-rhs;
lhs =y(30);
rhs =y(8);
residual(30)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(30, 30);

  %
  % Jacobian matrix
  %

  g1(1,1)=1-params(1)*(1-params(6));
  g1(1,2)=(-1)-(-(params(1)*(1-params(6))));
  g1(1,3)=(-T26);
  g1(1,4)=T26;
  g1(2,5)=1-params(4);
  g1(3,2)=1;
  g1(3,3)=1;
  g1(3,6)=params(9);
  g1(3,7)=(-1);
  g1(4,6)=(-1);
  g1(4,8)=1;
  g1(4,9)=(-1);
  g1(5,4)=(-1)-(-params(10));
  g1(5,10)=(-params(10));
  g1(5,22)=1;
  g1(6,11)=1;
  g1(6,12)=(-((1-params(11))/params(11)));
  g1(6,13)=(1-params(11))/params(11);
  g1(7,8)=(-(1-params(1)*params(11)));
  g1(7,11)=(-(params(1)*params(11)*params(12)));
  g1(7,12)=1-params(1)*params(11);
  g1(7,14)=(-(1-params(1)*params(11)));
  g1(8,8)=(-(1-params(1)*params(11)));
  g1(8,11)=(-(params(1)*params(11)*(params(12)-1)));
  g1(8,13)=1-params(1)*params(11);
  g1(9,7)=(-1);
  g1(9,9)=1;
  g1(9,14)=1;
  g1(9,15)=1;
  g1(10,1)=1;
  g1(10,2)=(-1);
  g1(10,4)=params(13);
  g1(10,10)=(-params(13));
  g1(11,2)=1;
  g1(11,15)=(-1);
  g1(11,16)=1;
  g1(11,20)=1;
  g1(12,2)=(-1);
  g1(12,15)=1;
  g1(12,16)=(-T109);
  g1(13,2)=1;
  g1(13,11)=(-1);
  g1(13,15)=(-1);
  g1(13,16)=T119;
  g1(13,23)=1;
  g1(14,11)=(-1);
  g1(14,24)=1;
  g1(15,2)=(-params(16));
  g1(15,3)=(-(params(3)*(1-params(14))/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,8)=1;
  g1(15,10)=(-((1-params(14))*params(18)/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,15)=params(16);
  g1(15,20)=(-(params(16)*params(16)*params(14)*params(19)*params(20)/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(15,25)=(-(params(16)*params(14)*params(19)*params(20)/(params(16)*params(14)*params(19)*params(20)+(1-params(14))*(params(3)+params(18)))));
  g1(16,1)=(-(T174*params(5)*params(7)));
  g1(16,2)=(-1)-(T158*(1-params(21)-(-params(21))-params(22))+T174*(-(params(5)*params(7))));
  g1(16,3)=(-(T174*params(3)*params(8)));
  g1(16,8)=1;
  g1(16,15)=1;
  g1(16,18)=(-(T158*(params(22)+params(21)-1-params(21))));
  g1(16,19)=(-(T158*params(22)));
  g1(16,22)=(-(T174*params(5)*params(7)));
  g1(17,19)=1;
  g1(17,26)=(-1);
  g1(18,1)=params(23);
  g1(18,4)=params(23);
  g1(18,5)=1;
  g1(18,18)=(-params(23));
  g1(18,19)=(-1);
  g1(18,21)=(-1);
  g1(18,27)=(-1);
  g1(19,21)=1;
  g1(20,19)=1-params(25);
  g1(21,23)=1;
  g1(22,4)=(-(1/params(6)*(1-(1-params(6)))));
  g1(22,22)=1;
  g1(23,9)=1-params(29);
  g1(24,27)=1-params(31);
  g1(25,24)=1-params(32);
  g1(26,25)=1-params(33);
  g1(27,26)=1-params(34);
  g1(28,18)=(-1);
  g1(28,28)=1;
  g1(29,2)=(-1);
  g1(29,29)=1;
  g1(30,8)=(-1);
  g1(30,30)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],30,900);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],30,27000);
end
end
end
end

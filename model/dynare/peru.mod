//Peru Model

var pd pc c d rl l w y a inv phi_h vn vd cm_h p_h ti tic q tau; //

varexo eps_a eps_v eps_phistar eps_ystar eps_rstar;

parameters  beta 
			psi 
			c_ss 
			eps 
			yd_ss 
			delta 
			pd_ss 
			pc_ss
			varphi
			alpha 
			omega 
			sigma_c 
			eta
			gamma
			ph_ss
			theta_h
			ti_ss
            inv_ss
			q_ss
			ystar_ss
			rl_ss
			r_ss 
			kappa 
			ttau
			gamma_r
			gamma_phi
			gamma_y
			gamma_b
            rho1a
            rho1v
            rho1p
            rho1y
            rho1r; 

beta= 0.99;
psi= 0.3;
c_ss= ;
eps= 0.75;
yd_ss= ;
delta= ;
pd_ss= ;
pc_ss= ;
varphi = 2;
alpha = ;
omega = ;
sigma_c = ;
eta = ;
gamma = ;
p_h_ss = ;
theta_h = ;
ti_ss = ;
inv_ss = ;
q_ss = ;
ystar_ss = ;
rl_ss = ;
r_ss = ;
kappa = ;
ttau = ;
gamma_r =  ;
gamma_phi = ;
gamma_y = ;
gamma_b = ; 
rho1a = ;
rho1v = ;
rho1p = ;
rho1y = ;
rho1r = ;



model(linear);
//1 consumption durable
pd-pc=(1-psi)/psi*(c_ss-eps*c_ss)/(yd_ss/delta)*pd_ss/pc_ss*((c-eps*c(-1))/(1-eps)-d)+beta*(1-delta)*(((c-eps*c(-1))-(c(+1)-eps*c))/(1-eps)+pd(+1)-pc(+1));

//2 consumption non durable
eps*(c-c(-1))= (c(+1)-c) -(1-eps)*(rl-(pc(+1)-pc));

//3 labor supply
varphi*l = w - pc ((c-eps*c(-1))/(1-eps)); 

//4 non-durable good prod. function
y = a + l;

//5 durable good prod. function
yd - d(-1) = alpha*(inv-d(-1));

//6 Home inflation
phi_h = (1-omega) / omega *(vn-vd);

//7
vn = (1-beta*omega)*(cm_h+y) +  beta*omega*(sigma_c*phi_h(+1) + vn(t+1));

//8
vd = (1-beta*omega)*y + beta*omega*((sigma_c-1)*phi_h(+1) + vd(t+1));

//9 marginal cost
cm_h = w - p_h - a;

//10
pd - pc = eta*(inv-d(-1));

//11 
q = -ti + p_h - pc;

//12
p_h - pc = gamma*(pc_ss/p_h_ss)^(theta_h-1)*ti;

//13
phic = phi_h + p_h(-1) - pc(-1) - gamma*(pc_ss*ti_ss/p_h_ss)^(theta_h-1)*ti;

//14
tic = ti(-1) + phi_h - phistar_f - e + e(-1);

//15 aggregate output
y= -theta_h*(p_h-pc) + ((1-gamma)*(c_ss*c+inv_ss*inv)+(gamma*q_ss*ystar_ss)*(theta_h*q + ystar))/((1-gamma)*(c_ss+inv_ss)+ gamma*q_ss*ystar_ss);

//16 bop 
p_h-pc+y = (b_ss/(p_h_ss*y_ss))*(rl_ss*(b(-1)+rl(-1)-pc) - b + pc - rl_ss*(b+rl-pc) + r_ss*(r+b-pc))  + (1/(p_h_ss*y_ss))*(pc_ss*c_ss*c + pd_ss*yd_ss*(pd-pc+yd));

//17 uip
r = rstar + e(+1) - e;
 
//18 lending interest rate
rl = r + v + kappa*(b-pd-d) + tau;

//19 macroprudential tool
tau = (b(-1) - b(-2))*ttau;

//20 taylor rule
r = gamma_r*r(-1) + (1-gamma_r)*( gamma_phi*(pc(-1)-pc(-2))  + gamma_y*(y(-1)-y(-2))  + gamma_b*(b(-1)-b(-2)));

//21 inflation
phic = pc - pc(-1);

//22
yd = 1/delta*(d- (1-delta)*d(-1));


//exogenous process
//23
a = rho1a*a(-1) + rho2a*a(-2)+ eps_a;

//24
v = rho1v*v(-1) + eps_v;

//25
phistar_f = rho1p*phistar_f(-1) + eps_phistar;

//26
ystar = rho1y*ystar(-1) + rho2y*ystar(-2) + eps_ystar;

//27
rstar = rho1r*rstar(-1) + rho2r*rstar(-2) + eps_rstar;

end;

shocks;
var eps_a; stderr exp(sigma_a);
var eps_v; stderr exp(sigma_v);
var eps_p; stderr exp(sigma_p);
var eps_y; stderr exp(sigma_y);
var eps_r; stderr exp(sigma_r);
end;


check;
steady;


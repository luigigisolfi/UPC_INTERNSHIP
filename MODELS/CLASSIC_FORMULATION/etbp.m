function xdot = etbp(~, x)
%------------------------------------------------------------------------
% Elliptical Three Body Problem vectorfield.
% Large Mass, 1-mu (Earth)    to the right of the origin at (mu, 0, 0)
% Small Mass, mu   (Moon) to the left at (mu-1, 0, 0).
%
%                      L5
%
% L2 -- Moon-- L1 ----------------- Earth --------------- L3
%
%                      L4
%
% Input variables: t (time), x (3D state, pos+vel) and mu (global mass param.)
% Output: f (vectorfield)
%-----------------------------------------------------------------------

global mu; 
global avg_e;
global nu

e= avg_e;

rtbp_pos_body_1 =  [mu;0;0];
rtbp_pos_body_2 = [(mu-1);0;0];

x_s1 = x(1:3)-rtbp_pos_body_1;
x_s2 = x(1:3)-rtbp_pos_body_2;

rho_s13 = norm(x_s1)^3;
rho_s23 = norm(x_s2)^3;

Ox = x(1) - ((1-mu)*x_s1(1)/rho_s13 + mu*x_s2(1)/rho_s23);
Oy = x(2) - ((1-mu)*x_s1(2)/rho_s13 + mu*x_s2(2)/rho_s23);
Oz =      - ((1-mu)*x_s1(3)/rho_s13 + mu*x_s2(3)/rho_s23);

Delta_Omega =[Ox;Oy;Oz];
z_versor = [0,0,1];

% the following is the expression as found in Beom Park et al 
%(paper: Assessment of dynamical models for transitioning from the
%Circular Restricted Three-Body Problem to an ephemeris
%model with applications, equation 12)
% Note you also need the 7th equation in this case, in order to account for
% the true anomaly evolution.

xdot = zeros(7,1);
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);
xdot(4:6) = - e*sin(nu)/(2*sqrt(1+e*cos(nu)))*xdot(1:3) - (2*cross(z_versor,xdot(1:3))*sqrt(1+e*cos(nu)))' - (e*cos(nu)*x(3)*z_versor)' + Delta_Omega;
xdot(7) = sqrt(1+avg_e*cos(nu)); %added this 7th element cause the true anomaly is also changing. this is, therefore, nu_dot.




end

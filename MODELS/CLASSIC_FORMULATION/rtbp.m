function xdot = rtbp(~, x)
%------------------------------------------------------------------------
% Spatial Circular RTBP vectorfield.
% Large Mass, 1-mu (Earth)    to the right of the origin at (mu, 0, 0)
% Small Mass, mu   (Moon) to the left at (mu-1, 0, 0).
%
%                      L5
%
% L2 -- Moon-- L1 ----------------- Earth --------------- L3
%
%                      L4
%
% Input variables: t (time), x (3D state, pos+vel) and mu (mass param.)
% Output: f (vectorfield)
%-----------------------------------------------------------------------
global mu;

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

xdot = zeros(6,1);

%this is equation 10 of Beom Park et al
% (Assessment of dynamical models for transitioning from the
%Circular Restricted Three-Body Problem to an ephemeris
%model with applications)
xdot(1:3) = x(4:6);
xdot(4:6) = -2*cross(z_versor, xdot(1:3))' + Delta_Omega;

end

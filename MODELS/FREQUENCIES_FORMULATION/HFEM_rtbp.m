function xdot = HFEM_rtbp(~,x)

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

%------------------------------------------------------------------------
% This is the same as rtbp.m vector field, [INSTRUCTIVE TO CHECK THAT THEY INDEED GIVE SIMILAR/EQUAL results!]
% only it is written in the framework used in Gomez et al 
% (Solar System with Selected
% Frequencies paper) or Beom Park et al (Assessment of dynamical models for transitioning from the
%Circular Restricted Three-Body Problem to an ephemeris
%model with applications), namely using the coefficients b1-13, 
% so as to be coherent with the HFEM vector field formulation. 
% Clearly, all coefficients are constant and do not depend on the true
% anomaly, nor eccentricity (since we're in the circular - not elliptical - case)
%-----------------------------------------------------------------------

global mu

b1=0;
b2=0;
b3=0;
b4=0;
b5 = 2;
b6=0;
b7 = 1;
b8=0;
b9=0;
b10 = 1;
b11=0;
b12=0;
b13 = 1;

xdot = zeros(6,1);

first = [b1;b2;b3];
second = [b4,b5,0;-b5,b4,b6;0,-b6,b4];
third = [b7,b8,b9;-b8,b10,b11;b9,-b11,b12];

rtbp_pos_body_1 =  [mu;0;0];
rtbp_pos_body_2 = [(mu-1);0;0];

x_s1 = x(1:3)-rtbp_pos_body_1;
x_s2 = x(1:3)-rtbp_pos_body_2;

rho_s13 = norm(x_s1)^3;
rho_s23 = norm(x_s2)^3;

synodic_acc_primaries = -(1-mu)*x_s1/rho_s13 - (mu)*x_s2/rho_s23; %primaries contribution to synodic acceleration

Delta_Omega = synodic_acc_primaries;

xdot(1:3) = x(4:6);
xdot(4:6) = first + second*x(4:6) + third*x(1:3) + b13*Delta_Omega;


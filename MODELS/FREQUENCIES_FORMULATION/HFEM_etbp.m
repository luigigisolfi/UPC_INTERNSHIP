function xdot = HFEM_etbp(~,x)
%------------------------------------------------------------------------
% Spatial Elliptical ETBP vectorfield.
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
% This is the same as etbp.m vector field, [INSTRUCTIVE TO CHECK THAT THEY INDEED GIVE SIMILAR/EQUAL results!]
% only it is written in the framework used in Gomez et al 
% (Solar System with Selected
% Frequencies paper) or Beom Park et al (Assessment of dynamical models for transitioning from the
%Circular Restricted Three-Body Problem to an ephemeris
%model with applications), namely using the coefficients b1-13, 
% so as to be coherent with the HFEM vector field formulation. 
% In this case, the coefficients depend on true anomaly and eccentricity of
% the primaries
%-----------------------------------------------------------------------
global mu
global avg_e
global nu


if ~exist('avg_e', 'var') || isempty(avg_e)
    error('Impossible to run: HFEM_etbp: Reason: global eccentricity variable (avg_e) does not exist or is empty.');
end

if ~exist('nu', 'var') || isempty(nu)
    error('Impossible to run: HFEM_etbp: Reason: global true anomaly variable (nu) does not exist or is empty.');
end

b1=0;
b2=0;
b3=0;
b4=-avg_e*sin(nu)/2*(sqrt(1+avg_e*cos(nu)));
b5 = 2*(sqrt(1+avg_e*cos(nu)));
b6=0;
b7 = 1;
b8=0;
b9=0;
b10 = 1;
b11=0;
b12=-avg_e*cos(nu);
b13 = 1;

xdot = zeros(7,1);

first = [b1;b2;b3];
second = [b4,b5,0;-b5,b4,b6;0,-b6,b4];
third = [b7,b8,b9;-b8,b10,b11;b9,-b11,b12];

rtbp_pos_body_1 = [(mu);0;0];
rtbp_pos_body_2 = [(mu-1);0;0];

x_s1 = x(1:3)-rtbp_pos_body_1;
x_s2 = x(1:3)-rtbp_pos_body_2;

rho_s13 = norm(x_s1)^3;
rho_s23 = norm(x_s2)^3;

synodic_acc_primaries = -(1-mu)*x_s1/rho_s13 - (mu)*x_s2/rho_s23; %primaries contribution to synodic acceleration

Delta_Omega = synodic_acc_primaries;

xdot(1:3) = x(4:6);
xdot(4:6) = first + second*x(4:6) + third*x(1:3) + b13*Delta_Omega;
xdot(7) = sqrt(1+avg_e*cos(nu)); %added this 7th element cause the true anomaly is also changing. this is, therefore, nu_dot.



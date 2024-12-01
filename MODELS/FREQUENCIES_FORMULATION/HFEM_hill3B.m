function xdot = HFEM_rtbp(t,x)

%------------------------------------------------------------------------
% This is supposed to be the same as rtbp vector field, [PLEASE CHECK THAT THEY INDEED GIVE SIMILAR/EQUAL results!]
% only it is written in the framework used in Gomez et al 
% (Solar System with Selected
% Frequencies pape), namely using the coefficients b1-13, 
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
b7 = 3; %coefficient for Hill
b8=0;
b9=0;
b10 = 0;
b11=0;
b12=-1; %coefficient for Hill
b13 = 1;

xdot = zeros(6,1);

first = [b1;b2;b3];
second = [b4,b5,0;-b5,b4,b6;0,-b6,b4];
third = [b7,b8,b9;-b8,b10,b11;b9,-b11,b12];

% Adjust position for Moon-centered frame
x(1) = x(1) - mu + 1;  % Shift to Moon-centered coordinates

% Compute gravitational acceleration terms
x_s2 = x(1:3);  % Position relative to the Moon (now the origin)
rho_s23 = norm(x_s2)^3;

synodic_acc_primaries = -(mu) * x_s2 / rho_s23;  % Gravitational acceleration

Delta_Omega = synodic_acc_primaries;

%vector field as derived in:
%https://core.ac.uk/download/pdf/20310212.pdf, section 3.9.1
%but using the b1-13 coefficient formulation
xdot(1:3) = x(4:6);
xdot(4:6) = first + second*x(4:6) + third*x(1:3) + b13*Delta_Omega;


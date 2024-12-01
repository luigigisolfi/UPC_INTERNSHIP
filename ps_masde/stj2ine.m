function  ene = stj2ine(x)
%-------------------------------------------------------------------------
% Energy of a state of the vectorfield stj2in (2BP, J2 perturbed)
% Input: x (3D state, pos+vel) 
%        Inside function: mu (mass param.) re (Earth radius) 
%             and xj2 (J2 value) all them in coherent units.
% Output: ene, the energy.
%-------------------------------------------------------------------------
% In canonical-like Earth units:
  mu=1.e0; re=1.e0; xj2=0.00108263e0;
%-------------------------------
aux=mu*re*re*xj2;
r2=x(1)*x(1)+x(2)*x(2);
z2=x(3)*x(3);
p2=1.e0/(r2+z2);
p1m=sqrt(p2);
epot=p1m*(aux*p2*(1.5*z2*p2-0.5)-mu);
ene=0.5*(x(4)*x(4)+x(5)*x(5)+x(6)*x(6))+epot;
end

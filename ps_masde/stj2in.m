function f = stj2in(t, x)
%-----------------------------------------------------------------------
% Vectorfield of the inertial Two Body Problem including the J2
% perturbation. (The main problem in artificial satellite theory).
%
% Input: t (time), x (3D state, pos+vel) 
%        Inside function: mu (mass param.) re (Earth radius) 
%             and xj2 (J2 value) all them in coherent units.
% Output: f (vectorfield)
%-----------------------------------------------------------------------
% In canonical-like Earth units:
  mu=1.e0; re=1.e0; xj2=0.00108263e0;
%-------------------------------
cta=1.5*mu*re*re*xj2;
r2=x(1)*x(1)+x(2)*x(2);
z2=x(3)*x(3);
p2=1.e0/(r2+z2);
p3m=p2*sqrt(p2);
aux=p3m*(mu+cta*p2*(1-5*z2*p2));
xdot=zeros(6,1);
xdot(1)=x(4);
xdot(2)=x(5);
xdot(3)=x(6);
xdot(4)=-aux*x(1);
xdot(5)=-aux*x(2);
xdot(6)=-x(3)*p3m*(mu+cta*p2*(1-p2*(3*z2-2*r2)));

f = xdot;

end

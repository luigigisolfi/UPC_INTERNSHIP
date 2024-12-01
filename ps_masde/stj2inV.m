function dy = stj2inV(t, x)
%-----------------------------------------------------------------------
% Vectorfield of the inertial Two Body Problem including the J2 
% perturbation (this is, the main problem in artificial satellite theory),
% INCLUDING VARIATIONAL equations.
%
% Input: t (time), x (3D state, pos+vel) 
%        Inside function: mu (mass param.) re (Earth radius) 
%             and xj2 (J2 value) all them in coherent units.
%
% Output: dy (vectorfield + variational equations in columns)
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
xdot=zeros(42,1);
xdot(1)=x(4);
xdot(2)=x(5);
xdot(3)=x(6);
xdot(4)=-aux*x(1);
xdot(5)=-aux*x(2);
xdot(6)=-x(3)*p3m*(mu+cta*p2*(1-p2*(3*z2-2*r2)));
o11=p3m*(p2*(3*mu*x(1)*x(1)-cta+p2*(5*cta*(x(1)*x(1)+z2)...
      -35*cta*x(1)*x(1)*z2*p2))-mu);
o22=p3m*(p2*(3*mu*x(2)*x(2)-cta+p2*(5*cta*(x(2)*x(2)+z2)...
      -35*cta*x(2)*x(2)*z2*p2))-mu);
aux1=21*z2-14*r2;
o33=p3m*(p2*(3*mu*z2-cta+cta*p2*(2*(7*z2-r2)-z2*p2*aux1))-mu);
o12=x(1)*x(2)*p3m*p2*(3*mu+p2*(5*cta-35*cta*z2*p2));
aux2=p3m*p2*(3*mu+p2*cta*(1-p2*aux1));
o13=x(1)*x(3)*aux2;
o23=x(2)*x(3)*aux2;
for j=7:6:37
  xdot(j)=x(j+3);
  xdot(j+1)=x(j+4);
  xdot(j+2)=x(j+5);
  xdot(j+3)=o11*x(j)+o12*x(j+1)+o13*x(j+2);
  xdot(j+4)=o12*x(j)+o22*x(j+1)+o23*x(j+2);
  xdot(j+5)=o13*x(j)+o23*x(j+1)+o33*x(j+2);
end

dy = xdot;

end

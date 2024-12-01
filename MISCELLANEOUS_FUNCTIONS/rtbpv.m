function xdot = rtbpv(t, x)
%-------------------------------------------------------------
% Spatial Circular BTBP vectorfield with variational equations.
% Large Mass, 1-mu (Earth) to the right of the origin at (mu, 0, 0)
% Small Mass, mu   (Moon) to the left at (mu-1, 0, 0).
%
%                L5
% L2- Moon -L1 ------------ Earth ---------- L3
%                L4
%
% Input variables: t (time), x (3D state, pos+vel+ini STM) and 
% mu (mass in global param.)
% Output: xdot (vectorfield)
%-------------------------------------------------------------
global mu;
xdot = zeros(42,1);

y1=x(1)-mu;
y12=y1*y1;
y2=x(2)*x(2);
y3=x(3)*x(3);
r1=y12+y2+y3;
r1a=r1*sqrt(r1);
r15=r1a*r1;
r1=r1a;
r2=(y1+1.e0)*(y1+1.e0)+y2+y3;
r2a=r2*sqrt(r2);
r25=r2a*r2;
r2=r2a;
p1=(1.e0-mu)/r1;
p2=mu/r2;
q=-1.e0*(p1+p2);

xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);
xdot(4) = 2*x(5) + x(1)-y1*p1-(y1+1.e0)*p2;
xdot(5) =-2*x(4) + x(2)*(1.e0+q);
xdot(6) = x(3)*q;

% Variational equations
rr1=3.e0*(1.e0-mu)/r15;
rr2=3.e0*mu/r25;
qp=rr1*y1+rr2*(y1+1.e0);
qpp=(rr1+rr2)*x(2);
o11=1.e0+q+rr1*y12+rr2*(y1+1.e0)*(y1+1.e0);
o12=qp*x(2);
o13=qp*x(3);
o22=1.e0+q+(rr1+rr2)*y2;
o23=qpp*x(3);
o33=q+(rr1+rr2)*y3;
for j6 = 6:6:36
  xdot(j6+1:j6+3)=x(j6+4:j6+6);
  xdot(j6+4)=o11*x(j6+1)+o12*x(j6+2)+o13*x(j6+3)+2.e0*x(j6+5);
  xdot(j6+5)=o12*x(j6+1)+o22*x(j6+2)+o23*x(j6+3)-2.e0*x(j6+4);
  xdot(j6+6)=o13*x(j6+1)+o23*x(j6+2)+o33*x(j6+3);
end
end

   function cj = cjrtbp(x,mu)
%----------------------------------------------------------
% Computes the Jacobi constant of a spatial RTBP state x
%----------------------------------------------------------
 x12=x(1)*x(1); x22=x(2)*x(2); x32=x(3)*x(3);
 r1=sqrt((x(1)-mu+1)^2+x22+x32);
 r2=sqrt((x(1)-mu)^2+x22+x32);
 Ome=0.5*(x12+x22)+mu/r1+(1-mu)/r2+0.5*mu*(1-mu);
 cj=2.0*Ome-x(4)*x(4)-x(5)*x(5)-x(6)*x(6);
end

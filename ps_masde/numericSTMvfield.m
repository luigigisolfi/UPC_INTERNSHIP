function STM = numericSTMvfield(ti,tf,xi,eps,vfield,hmin,hmax,tol)
%------------------------------------------------------------------------
% By means of numerically differentiation it computes the  State Transition 
% Matrix from ti to tf of an initial condition xi in the vectorfield
% vfield.
% n is the dimension of the vectorfield and eps is the step to be used for 
% the numerical differentiation. hmin, hmax, tol are resp. min, max and 
% truncation error for the numerical propagation.
% The STM returned is a n*n matrix.
%------------------------------------------------------------------------
STM=zeros(6,6);
for k=(1:6)
  xa=xi;
  xa(k)=xi(k)+eps;
  xfp=propTITF_vfield(ti,xa,tf,vfield,hmin,hmax,tol);
  xa(k)=xi(k)-eps;
  xfm=propTITF_vfield(ti,xa,tf,vfield,hmin,hmax,tol);
  STM(:,k)=(xfp-xfm)/(2*eps);
end
end


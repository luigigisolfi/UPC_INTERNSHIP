% ------------------------------------------------------------------------
% Checks the propagated State Transition Matrix against the one obtained
% by numerical differentiation. (Vectorfield checker)
% ------------------------------------------------------------------------
clear all; close all; clc;
%--------------------------------------------------------------------------
 ti=0;
 xi=[-1.2499884859E-00  .0000000000E+00  1.0000000000E+00 ...
     .0000000000E+00    .9163258721E-01 .0000000000E+00];
 tf=1.e0;
 eps=1.e-5;
 hmin=1.e-4;
 hmax=1.e0;
 tol=1.e-10;
%--------------------------------------------------------------------------
xiv=zeros(1,42);
xiv(1:6)=xi;
for i= 1:6 , xiv(7*i)=1.e0; end  % variational matrix initalized with Id.
xf=propTITF_vfield(ti,xi,tf,@stj2in,hmin,hmax,tol);
xfv=propTITF_vfield(ti,xiv,tf,@stj2inV,hmin,hmax,tol);
stm=numericSTMvfield(ti,tf,xi,eps,@stj2in,hmin,hmax,tol);
eni=stj2ine(xi);
enf=stj2ine(xf);
fprintf('Energy. initial, final and diff: %e %e %5.2e\n',eni,enf,enf-eni);
for i = 0:5,  stmv(6*i+1 : 6*i+6)=stm(:,i+1); end
for i= 1:36 
fprintf('diff in STM: %2d %14.6e %14.6e %10.2e\n',i+6,xfv(i+6),stmv(i),xfv(i+6)-stmv(i));
end

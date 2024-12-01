function [tk1,xk1,hk1,err] = rk45f(tk,xk,hk,hmin,hmax,tol,vfield)
%-------------------------------------------------------------------------
% Runge-Kutta Fehlberg orders 4-5.
% Input: tk (time), xk (state), hk (suggested time step forward (+) or
%        backward (-), hmin (min allowed step), hmax (max allowed step),
%        tol (truncation error), vfield (vectorfield).
% Output: tk1, xk1: time and state forward or backward in the propagation
%         hk1: suggested time step for the next call
%         err: truncation error found in the propagation from tk to tk1
%
% Note. It accepts a call with abs(hk)<hmin. In this case it proceeds to
%       integrate with step h regardless of error.
%-------------------------------------------------------------------------
[nr,nc]=size(xk);
if (nr<nc), xk=xk'; end % input xk can be either row or column vector
error=realmax; hs=hk;
if (abs(hk)>hmin), is=0; else, is=1; end
while (error > tol && is<2)
 if (is==1), is=2; end  % if is=1 this will be last iteration (min step)
 ha=hs;
 k1=ha*vfield(tk,xk);
 k2=ha*vfield(tk+0.25*ha,xk+0.25*k1);
 k3=ha*vfield(tk+0.375*ha,xk+0.09375*k1+0.28125*k2);
 k4=ha*vfield(tk+(12/13)*ha,xk+(1932/2197)*k1-(7200/2197)*k2+(7296/2197)*k3);
 k5=ha*vfield(tk+ha,xk+(439/216)*k1-8*k2+(3680/513)*k3-(845/4104)*k4);
 k6=ha*vfield(tk+ha/2,xk-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-0.275*k5);
 error=norm(k1/360-(128/4275)*k3-(2197/75240)*k4+k5/50+(2/55)*k6);
 hp=abs(ha/hmax);
 if error<=tol*hp*hp*hp*hp*hp
   hs=hmax*sign(ha);
 else
   hs=0.9*ha*(tol/error)^0.2;
   if (abs(hs)<=hmin), hs=hmin*sign(ha); is=is+1; end
 end  
end
err=error;
hk1=hs;
tk1=tk+ha;
xk1=xk+(16/135)*k1+(6656/12825)*k3+(28561/56430)*k4-(9/50)*k5+(2/55)*k6;
if (nr<nc), xk1=xk1'; end
end

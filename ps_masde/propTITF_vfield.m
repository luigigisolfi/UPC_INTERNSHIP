function xf = propTITF_vfield(ti,xi,tf,vfield,hmin,hmax,tol)
%--------------------------------------------------------------------------
% Propagates an inital condition xi at time ti up to time tf in
% the vectorfield vfield. hmin, hmax, tol are the parameters for rk45f. 
% NOTE: initial step h is set in the first line of the function.
%--------------------------------------------------------------------------
 h=1.e-2;
 xf=xi;
 if (abs(ti-tf) < 1.e-12), return; end
 if (tf < ti), h=-h; end
 t=ti;
 while ((t<tf & h>0) || (t>tf & h<0)) 
 [t,xf,h,err]=rk45f(t,xf,h,hmin,hmax,tol,vfield);
 end
 while (abs(t-tf) > 1.e-12)
 h=tf-t;
 [t,xf,h,err]=rk45f(t,xf,h,hmin,hmax,tol,vfield);
 end
end
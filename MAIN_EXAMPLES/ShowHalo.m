%--------------------------------------------------------------------------
% This is a script given to me by professor Masdemont
% showcasing some initial halo orbit conditions for different values of mu.
% I used it to validate my propagators: if the propagator is correct, 
% I should get a halo orbit. 
% Unfortunately, no initial conditions for mu (sun-earth) nor mu
% (earth-moon) are available in this script.
%--------------------------------------------------------------------------
clear all;
global mu;
 hmin=1.e-6; hmax=1.e0; tol=1.e-10; h=0.001;
 ntall=1; iorb=-1;
%------------------------------ SELECT ORBIT HERE with variable ncas 
ncas=3;
switch(ncas)
  case 1
   mu=0.1;
   xi=[-6.356279410857e-01 0.e0 2.80e-2   0.e0 2.450072408321e-01 0.e0];
   Period=2.449380491547e+00;  
  case 2
   mu=0.1;
   xi=[-6.426091998376e-01 0.e0 5.58e-2  0.e0 2.828055215743e-01 0.e0];
   Period=2.465684027325e+00;   
  case 3
   mu=0.2;
   xi=[-4.691774486730e-01 0.e0 5.30e-2   0.e0 3.086553889699e-01 0.e0];
   Period=2.345047489718e+00;
  case 4
   mu=0.2;    
   xi=[-4.797292756212e-01 0.e0 8.78e-2  0.0e0 3.717293897424e-01 0.e0];
   Period=2.375609424166e+00;
end

%-------------------------------  PROPAGATION
ti=0; tf=10*Period;
h=0.01; hmin=1.e-6; hmax=1; tol=1.e-12; iorb=1;
[tfa,xf]= propTITF(ti,xi,tf/5,@HFEM_etbp,hmin,hmax,tol,iorb);
% ------ plot orbit
plot3(xf(:,1),xf(:,2),xf(:,3))

%------- accuracy of propagation
cjf=cjrtbp(xf(end,:),mu); cji=cjrtbp(xi,mu);
fprintf('------------- Result of integration ------------------ \n');
fprintf('Value of mu: %f, Period: %19.12e\n',mu,Period);
fprintf('IC   (x,y,z): %19.12e %19.12e %19.12e \n',xi(1:3));
fprintf('  (xd,yd,zd): %19.12e %19.12e %19.12e \n',xi(4:6));
fprintf('CT Jacobi: %19.12e. Dif initial-final: %10.3e\n',cji,cjf-cji);
fprintf('Norm of diference states (p+v) initial-final OP: %e\n',norm(xi(1:6)-xf(1:6)));

% Combine data into a matrix
data = [tfa, xf(:,1), xf(:,2), xf(:,3), xf(:,4), xf(:,5), xf(:,6)];
% Open a file for writing
fileID = fopen('orbit_data.txt', 'w');

% Check if the file opened successfully
if fileID == -1
    error('Failed to open the file.');
end

% Write the data to the file
% The format string specifies the format for each column
formatSpec = '%f %f %f %f %f %f %f\n';
fprintf(fileID, formatSpec, data.');

% Close the file
fclose(fileID);

disp('Data successfully written to orbit_data.txt');









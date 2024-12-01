clear all; close all; clc;
%--------------------------------------------------------------------------
% This is a main example to compute the residual acceleration score for
% different models. 
% This code:
% 1) reads an orbit (output of perform_parallel_shooting function)
% 2) computes the vector field values along the orbit (HFEM_RTBP, HFEM_ETBP, HFEM_ETBP_RESIDUAL)
% 3) compute the difference of HFEM - HFEM_model for each of the 3 models
% 4) for each point along the orbit, it keeps track of which model is
% "closest" to HFEM. Circles: ETBP, Crosses: RTBP, Squares: ETBP_RESIDUAL
% 5) assigns a normalized score 
% 6) the points of the orbit are plotted together with Circles, Crosses or
% Squares, in order to give a visual representation of how "good" each
% model performs along the trajectory. 

%NOTE: 
% this was specifically done for a Earth-Moon Lissajous orbit, and the SUN
% is the only perturbation. Jupiter can be added easily in the list of
% bodies as well, but in that case you also need to tweak the commented lines:
% 158-159-160.

%NOTE:
%In this example, we compute - among others - the vector field
%HFEM_etbp_residual. This is basically HFEM - HFEM_etbp. 
%In other words, instead of computing HFEM(t,x) - HFEM_etbp(t,x) (as it is
%done in the etbp_rtbp_hilleb_lissajous_comparison_3d_orbit example), 
%we  compute directly: HFEM_etbp_residual(t,x) thus reducing the computational time. 
%Of course, this can be done with any other vector field, by just taking
%the differences of the b1-13 coefficients. For instance, 
% I strongly advice you to write the HFEM_rtbp_residual and
% HFEM_hill3b_residual vector fields, just as I did with
% HFEM_etbp_residual. 

%NOTE:
% This example might not give the exact same output as its "twin example", 
% since it was written just as a test to showcase the potential to evaluate
% differences of vector fields directly as it is done in
% HFEM_etbp_residual. 
% SO DONT WORRY IF YOU SEE SOME ERRORS IN OUTPUT :)
%--------------------------------------------------------------------------
global avg_e
global interpolators
global n_rtbp
global n_anomalistic
global PRIMARIES
global mu
global BODIES
global L
global eclipse_date_et
global nu


% --------------- LOAD KERNELS -------------
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------- %

%  --------------- SET THE MODEL -------------
FRAME = 'J2000';
OBSERVER = 'SUN';
%BODIES = [{'JUPITER BARYCENTER'}, {'MARS BARYCENTER'},{'MOON'}, {'SATURN BARYCENTER'}, {'MERCURY'}, {'VENUS'}, {'URANUS BARYCENTER'}, {'PLUTO BARYCENTER'}, {'NEPTUNE BARYCENTER'}];
BODIES = {'SUN'};
PRIMARIES = [{'EARTH'},{'MOON'}];
L = get_L(PRIMARIES, FRAME); % Computes the mean distance of two given primaries over a timespan of 50 years 
[mu, body_1, body_2] = get_mu(PRIMARIES); %Compute grav. parameter for the system
MODEL = '@etbp'; 
%-------- Compute mean anomaly (global variable, needed for time converison) ---------------%
%n_anomalistic = 2.639394888546285e-06 %computed with the function n_peri_peri, corresponding to the anomalistic month
n_anomalistic = get_n_peri_peri(PRIMARIES, FRAME); %computed over 50 years timespan
n_rtbp = get_n(PRIMARIES,FRAME); %computed with Kepler's 3rd Law
%------------------------------------------------------------------------------%

%-----Set LUNAR ECLIPSE of 21 JAN 2000 as starting point in time (taken from almanacs and verified with SPICE) -----%
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
eclipse_date_UTC= '21 JAN 2000 04:05:02';
eclipse_date_et = cspice_str2et(eclipse_date_UTC);
t_list = linspace(eclipse_date_et,eclipse_date_et + 100*pi/n_rtbp, 27*100);
%this will be needed later
[inertial_state_primaries, inertial_state_bodies, interpolators] = get_ephemeris(t_list, PRIMARIES, BODIES, FRAME, OBSERVER);

cspice_kclear()
%-------------------------------------------------------------------------------------------------------------%

% ------ Critical Part: Defining the ephemeris times (inertial, common for all models) ------ %
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
ti_utc = cspice_et2utc(t_list(1), 'C', 0);
tf_utc = cspice_et2utc(t_list(end), 'C', 0);
cspice_kclear()

fprintf('-----------------------------------------------------------\n')
fprintf('MAIN\nUTC Start Simulation Epoch: %s\n', ti_utc)
fprintf('UTC End Simulation Epoch: %s\n', tf_utc)
fprintf('-------------------------------------------------------------------------------------------------------------------\n')
% -----------------------------------------------------------------%

%--------- Retrieve INERTIAL state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE) ---------%
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
%[inertial_state_primaries, inertial_state_bodies, interpolators] = get_ephemeris(t_list, PRIMARIES, BODIES, FRAME, OBSERVER);
avg_e = get_avg_e(PRIMARIES, FRAME); %computes average eccentricity over the t_list timespan. Done by interpolating ephemeris
nu = get_nu_eclipse(PRIMARIES, eclipse_date_et);
% ------------------------ RTBP variables ---------------------------%
GM_1 = get_GM_body(body_1);
GM_2 = get_GM_body(body_2);
%---------------------------------------------------------------------%

% ------ Pick pick a test orbit (output of parallel shooting) --------- %
%Test Orbit
filename = './ORBITS/test_lissajous_for_comparison.txt';
[t_orb,x_orb] = read_orbit(filename);
%---------------------------------------------------------------------%

% ---- Initialize and compute acceleration values for contour plot ----- %
Z = zeros(size(x_orb(:,1)));
Z_rtbp = zeros(size(x_orb(:,1)));
Z_etbp = zeros(size(x_orb(:,1)));
Z_residual = zeros(size(x_orb(:,1)));

vector_field = zeros(3, length(x_orb(:,1)));

% Define a state vector for each point
for i = 1:1:length(x_orb(1:1000,1))
        i
        t_HFEM = t_orb(i)/n_rtbp;
        t_peri_peri = (t_HFEM - eclipse_date_et)*n_anomalistic;
        t_rtbp = (t_HFEM - eclipse_date_et)*n_rtbp;
        x = [x_orb(i,1), x_orb(i,2), x_orb(i,3), x_orb(i,4), x_orb(i,5), x_orb(i,6)].'; % Example state vector, modify as needed
        
        % Compute the vector fields
        f1_val = HFEM(t_HFEM*n_anomalistic, x);
        f2_val_full = HFEM_etbp(t_peri_peri, x); %translation of time needed: we want etbp to start from zero at t_HFEM. Note, however, the avg_nu_dot conversion
        f2_val = f2_val_full(1:6);
        f3_val = HFEM_rtbp(t_rtbp, x);        
        f4_val_full = HFEM_residual_etbp((t_HFEM)*n_anomalistic,x);
        f4_val = f4_val_full(1:6);

        %nu_dot_res = f4_val_full(7); %HFEM residual etbp gives nu as output, besides the vector state (6 + 1 elements)
        nu_dot = f2_val_full(7); %HFEM etbp b coefficients also gives the nu value. Here, we compute them to compare them

        Z_rtbp(i) =  abs(norm(f1_val(4:6) - f3_val(4:6))); %Difference of magnitudes of the acceleration vectors HFEM - RTBP
        Z_etbp(i) =  abs(norm(f1_val(4:6) - f2_val(4:6)));%Difference of magnitudes of the acceleration vectors HFEM - ETBP
        Z_residual(i) = abs(norm(f1_val(4:6) - f4_val(4:6)));%Difference of magnitudes of the acceleration vectors HFEM - ETBP RESIDUAL

        %Z_etbp(i) =  (norm(f4_val(4:6)));%Difference of magnitudes of the acceleration vectors HFEM - ETBP
        
        vector_field(:,i) = f1_val(4:6);
        nu = nu_dot*t_peri_peri + nu; %update the true anomaly assuming linearity (works quite well, same was done in Beom Park et al)

        %nu_res = nu_dot_res*t_peri_peri + nu;
end

 % ------- Initialize and compute index matrix, which will tell us if Z_rtbp >= Z_etbp or viceversa ----- %
index_matrix = zeros(size(Z_rtbp));
for i = 1:size(Z_rtbp, 1)
        if Z_rtbp(i) < Z_etbp(i) && Z_rtbp(i) < Z_residual(i)
            Z(i) = Z_rtbp(i);
            index_matrix(i) = 1;  % Indicates the value is from Z_rtbp
        elseif Z_etbp(i) < Z_rtbp(i) && Z_etbp(i) < Z_residual(i)
            Z(i) = Z_etbp(i);
            index_matrix(i) = 2;  % Indicates the value is from Z_etbp
        elseif Z_residual(i) < Z_rtbp(i) && Z_residual(i) < Z_etbp(i)
            Z(i) = Z_residual(i);
            index_matrix(i) = 3;
        else 
            Z(i) = Z_etbp(i);
            index_matrix(i) = 4;
        end
    if isnan(Z(i))
        fprintf('found a NaN at coordinates: %d, %d, %d.\n', x_orb(i,1),x_orb(i,2),x_orb(i,3))
        Z(i) = 0;
    end

end
% ------------------------------------------------------------------------------------------------%

% ----- Retrieve Additional Bodies Positions ---------------------------------------%
%inertial_pos_j =zeros(3,1);
%inertial_pos_j(1) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{1}, t_HFEM);
%inertial_pos_j(2) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{2}, t_HFEM);
%inertial_pos_j(3) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{3}, t_HFEM);
inertial_pos_m =zeros(3,1);
inertial_pos_m(1) =  ppval(interpolators.('SUN').spline{1}, t_HFEM);
inertial_pos_m(2) =  ppval(interpolators.('SUN').spline{2}, t_HFEM);
inertial_pos_m(3) =  ppval(interpolators.('SUN').spline{3}, t_HFEM);
%----------------------------------------------------------------------------------------%

% ----------------Evaluate the splines for the primaries
% ------------------------%
primary_name = PRIMARIES{1};
secondary_name = PRIMARIES{2};
primary_str = regexprep(primary_name, [{'\s+'}, {'-'}], '_');
secondary_str = regexprep(secondary_name, [{'\s+'}, {'-'}], '_');

% Initialize arrays to store interpolated position and velocity
rp = zeros(3, 1);
vp = zeros(3, 1);

rs = zeros(3, 1);
vs = zeros(3, 1);


for dim = 1:12
    interpolated_p = ppval(interpolators.(primary_str).spline{dim}, t_HFEM);
    interpolated_s = ppval(interpolators.(secondary_str).spline{dim}, t_HFEM);
    if dim <= 3
        rp(dim, :) = interpolated_p;
        rs(dim,:) = interpolated_s;
    elseif dim>=4 && dim<=6 
        vp(dim-3, :) = interpolated_p;
        vs(dim-3, :) = interpolated_s;
    elseif dim>=7 && dim<=9
        ap(dim-6, :) = interpolated_p;
        as(dim-6, :) = interpolated_s;       
    else
        oap(dim-9, :) = interpolated_p;
        oas(dim-9, :) = interpolated_s;    
    end
end
%----------------------------------------------------------------%

%retrieving primaries relative acceleration at time t
rs_rp = rs-rp;
vs_vp = vs-vp;
as_ap = as-ap;
oas_oap = oas-oap;

SEb_pos = rp + mu*rs_rp;
SEb_vel = vp + mu*vs_vp;
SEb_acc = ap + mu*as_ap;

b = SEb_pos;
b_dot = SEb_vel;
b_ddot = SEb_acc;
k = norm(rs_rp);

C = construct_C(rs_rp, vs_vp);

%rtbp_pos_j =  C\(inertial_pos_j-b)/k;
rtbp_pos_m =  C\(inertial_pos_m-b)/k;
% Plot the difference as a contour plot

figure;
hold on
[X,Y] = meshgrid(x_orb(:,1), x_orb(:,2));
% create corresponding Z values, assume z = 0 for locations with no z data
Z_grid = zeros(length(x_orb(:,1)),length(x_orb(:,2))) ;
for i = 1:length(x_orb(:,1))
    for j = 1:length(x_orb(:,2))
        if i==j % z data exist for only for x(n) y(n) location, n = 1,2,3...
        Z_grid(i,j) = Z(i);
        end
    end
end

contourf(X,Y,Z_grid)
title('Contour plot of interpolated data');
xlabel('x');
ylabel('y');
colorbar;

% Find the indices where index_matrix is 1
[row] = find(index_matrix == 1);

% Extract the corresponding X, Y, and Z values
crossX1 = x_orb(row,1);
crossX2 = x_orb(row,2);
% Overlay the crosses on the contour plot
scatter(crossX1, crossX2, 'xk','DisplayName', 'RTBP');

[row_] = find(index_matrix == 2);

% Extract the corresponding X, Y, and Z values
crossX1_ = x_orb(row_,1);
crossX2_ = x_orb(row_,2);
% Overlay the crosses on the contour plot
scatter(crossX1_, crossX2_,'o', 'DisplayName', 'ETBP');

[row__] = find(index_matrix == 3);
% Extract the corresponding X, Y, and Z values
crossX1__ = x_orb(row__,1);
crossX2__ = x_orb(row__,2);
% Overlay the crosses on the contour plot
scatter(crossX1__, crossX2__, 's', 'DisplayName', 'residual');

[row___] = find(index_matrix == 4);
% Extract the corresponding X, Y, and Z values
crossX1___ = x_orb(row___,1);
crossX2___ = x_orb(row___,2);
% Overlay the crosses on the contour plot
%scatter(crossX1__, crossX2__, 's', 'DisplayName', 'SAME');
legend('show')
return

%quiver3(mu-1, 0, 0, rtbp_pos_j(1)/(norm(rtbp_pos_j)*100), rtbp_pos_j(2)/(norm(rtbp_pos_j)*100), 0, 'r', 'LineWidth', 2, 'DisplayName','JUPITER BARYCENTER');
quiver3(mu-1, 0, 0, rtbp_pos_m(1)/(norm(rtbp_pos_m)*100), rtbp_pos_m(2)/(norm(rtbp_pos_m)*100), 0, 'y', 'LineWidth', 2, 'DisplayName','SUN DIRECTION');

for g = 1:(length(X))
    for d = 1:length(X)
        quiver3(X(g,d), Y(g,d),0, vector_field(1,g)*0.0001/(norm(vector_field(:,g))),vector_field(2,g)*0.0001/norm(vector_field(:,g)),0, 'HandleVisibility','off', 'LineWidth',1, 'Color','red')
    end
end
return


scatter(mu-1,0, 'filled', 'DisplayName','MOON')
th = 0:pi/50:2*pi;
x_hill = sqrt((mu)/3) * cos(th) + mu-1;
y_hill = sqrt((mu)/3) * sin(th) + 0;
plot(x_hill,y_hill, 'LineWidth',2)

if strcmp(MODEL, '@etbp')
    title_str = sprintf('a - difference between HFEM and ETBP (t = %s)', t_utc(1:11));
else
    title_str = sprintf('a - difference between HFEM and RTBP (t = %s)', t_utc(1:11));
end

title(title_str)
xlabel('x (adimensional)');
ylabel('y (adimensional)');
legend('show')
% Save the figure as a PNG
fileout = sprintf('%s.png', t_utc(1:11)); % 'plot_1.png', 'plot_2.png', ...
saveas(gcf, fileout); % Save as PNG

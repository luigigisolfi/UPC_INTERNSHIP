clear all; close all; clc;

global avg_e
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
%BODIES = [{'JUPITER BARYCENTER'}, {'MARS BARYCENTER'},{'SUN'}, {'SATURN BARYCENTER'}, {'MERCURY'}, {'VENUS'}, {'URANUS BARYCENTER'}, {'PLUTO BARYCENTER'}, {'NEPTUNE BARYCENTER'}];
BODIES = [{'MOON'}, {'JUPITER BARYCENTER'}];
PRIMARIES = [{'EARTH'}, {'MOON'}];
L = get_L(PRIMARIES, FRAME); % Computes the mean distance of two given primaries over a timespan of 50 years 
%L = 384601.25606767; % this is the value from Gomez et al (earth-moon)
[mu, body_1, body_2] = get_mu(PRIMARIES); %Compute grav. parameter for the system
%-------- Compute mean anomaly (global variable, needed for time converison) ---------------%
n_anomalistic = get_n_peri_peri(PRIMARIES, FRAME);
n_rtbp = get_n(PRIMARIES,FRAME); %computed with Kepler's 3rd Law]
%------------------------------------------------------------------------------%

%-----Set LUNAR ECLIPSE of 21 JAN 2000 as starting point in time (taken from almanacs and verified with a SPICE plot) -----%
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
eclipse_date_UTC= '21 JAN 2000 04:05:02';
eclipse_date_et = cspice_str2et(eclipse_date_UTC);
avg_e = get_avg_e(PRIMARIES, FRAME); %computes average eccentricity over the t_list timespan. Done by interpolating ephemeris
nu = get_nu_eclipse(PRIMARIES, eclipse_date_et);
cspice_kclear()
%-------------------------------------------------------------------------------------------------------------%

% ------ Critical Part: Defining the ephemeris times (inertial, common for all models) ------ %
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
ti_utc = cspice_et2utc(eclipse_date_et, 'C',0);
tf_utc = cspice_et2utc(eclipse_date_et + 100*pi/n_rtbp, 'C', 0);
cspice_kclear()

fprintf('-----------------------------------------------------------\n')
fprintf('MAIN\nUTC Start Simulation Epoch: %s\n', ti_utc)
fprintf('UTC End Simulation Epoch: %s\n', tf_utc)
fprintf('-------------------------------------------------------------------------------------------------------------------\n')

%--------- Retrieve INERTIAL state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE) ---------%
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels

% ------------------------ RTBP variables ---------------------------%

GM_1 = get_GM_body(body_1);
GM_2 = get_GM_body(body_2);

% Initial conditions for satellite in Earth-Moon case
%initial_state_etbp = [mu - 1 + (mu/3)^(1/3); 0.01; 0.02; 0; -0.1; 0.02; nu];
%initial_state_rtbp = [mu - 1 + (mu/3)^(1/3); 0.01; 0.02; 0; -0.1; 0.02];

% Initial conditions for satellite in Sun-Earth case (close to L1 Earth)
%initial_state_etbp = [mu - 1 + (mu/3)^(1/3); 0.0; 0.0; 0; -0.01; 0.0; nu];
%initial_state_rtbp = [mu - 1 + (mu/3)^(1/3); 0.0; 0.0; 0; -0.01; 0.0];

% Initial conditions for satellite in Sun-Earth case (very close to Earth)
initial_state_etbp = [mu - 1 + (mu/60)^(1/3); 0.0; 0.0; 0; -0.01; 0.0; nu];
initial_state_rtbp = [mu - 1 + (mu/60)^(1/3); 0.0; 0.0; 0; -0.01; 0.0];

tspan = [0,3.14/5];
tspan_hfem = (tspan - eclipse_date_et*n_anomalistic);
tspan_peri_peri = (tspan/n_rtbp - eclipse_date_et)*n_anomalistic;

% Integrate the RTBP equations of motion
%[t_rtbp, state_rtbp] = ode78(@HFEM_rtbp, tspan_rtbp, initial_state_rtbp);
[t_hfem, state_hfem] = ode78(@HFEM, tspan_hfem, initial_state_rtbp);

% Integrate the ETBP equations of motion
[t_etbp, state_etbp] = ode78(@HFEM_etbp, tspan_peri_peri, initial_state_etbp);

% Plot the orbits in a 3D view for comparison
figure;
%plot3(state_rtbp(:,1), state_rtbp(:,2), state_rtbp(:,3), 'b', 'LineWidth', 1.5); % RTBP trajectory
hold on;
plot3(state_hfem(:,1), state_hfem(:,2), state_hfem(:,3), 'r', 'LineWidth', 1.5); % HFEM trajectory
plot3(state_etbp(:,1), state_etbp(:,2), state_etbp(:,3), 'b--', 'LineWidth', 1.5); % ETBP trajectory
%plot3(mu, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Primary 1
plot3(mu-1, 0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g'); % Primary 2
%earth_moon_sphere('moon_dist', 'both')
title('Comparison of Satellite Orbits: HFEM vs. ETBP');
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('HFEM Orbit', 'ETBP Orbit', 'Earth');
grid on;
axis equal;

% Plot the differences in position over time
position_diff = vecnorm(state_rtbp(:,1:3) - state_hfem_interp, 2, 2);
figure;
plot(t_rtbp, position_diff, 'k', 'LineWidth', 1.5);
title('Position Difference Between RTBP and ETBP Over Time');
xlabel('Time');
ylabel('Position Difference');
grid on;

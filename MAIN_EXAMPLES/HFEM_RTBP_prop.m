%------------------------------------------------------------------------
% Large Mass, 1-mu (Earth)    to the right of the origin at (mu, 0, 0)
% Small Mass, mu   (Moon) to the left at (mu-1, 0, 0).
%
%                      L5
%
% L2 -- Moon-- L1 ----------------- Earth --------------- L3
%
%                      L4
%
%-----------------------------------------------------------------------

% Main script to propagate satellite orbit using HFEM_rtbp
global mu
% Set global variables
mu = 0.012150582; % Example value for Earth-Moon system mass ratio

% Define initial state of the satellite (example values)
% [x0, y0, z0, vx0, vy0, vz0]
initial_state = [mu - 1 + (mu/3)^(1/3); 0.01; 0.02; 0; -0.1; 0.02]; % Position and velocity % For RTBP (6 elements)

% Define time span for the propagation
tspan = [0, 3]; % Start time and end time in nondimensional units

% Integrate the equations of motion
[t, state] = ode78(@HFEM_rtbp, tspan, initial_state);

% Plot the orbit
figure;
plot3(state(:,1), state(:,2), state(:,3), 'b'); % 3D plot of position
hold on;
earth_moon_sphere('moon_dist','both') %this is just for fun, but you can also just plot the primary at mu and secondary at mu-1
title('Satellite Orbit Propagation in the Synodic Frame');
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('Satellite Orbit RTBP', 'Earth', 'Moon');
grid on;
axis equal;

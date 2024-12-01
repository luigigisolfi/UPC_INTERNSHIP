% Main script to propagate orbit using HFEM_etbp_b_coefficients

% Set global parameters
global mu avg_e nu

% Define system-specific parameters (e.g., Earth-Moon system)
mu = 0.012150582; % Example: Earth-Moon mass ratio
avg_e = 0.0549; % Example: average eccentricity of the Moon’s orbit around Earth
nu = 0; % Initial true anomaly

% Define initial state of the satellite
% [x0, y0, z0, vx0, vy0, vz0, nu]
initial_state = [mu - 1 + (mu/3)^(1/3); 0.01; 0.02; 0; -0.1; 0.02; nu]; % Position and velocity % Example position, velocity, and initial anomaly
% Define time span for propagation
tspan = [0, 3]; % Time span in nondimensional units

% Integrate the equations of motion
[t, state] = ode78(@HFEM_etbp_b_coefficients, tspan, initial_state);

% Plot the orbit
figure;
plot3(state(:,1), state(:,2), state(:,3), 'b'); % 3D plot of position
hold on;
earth_moon_sphere('moon_dist', 'both')
title('Satellite Orbit Propagation with Elliptical Primaries in the Synodic Frame');
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('Satellite Orbit ETBP', 'Earth', 'Moon');
grid on;
axis equal;
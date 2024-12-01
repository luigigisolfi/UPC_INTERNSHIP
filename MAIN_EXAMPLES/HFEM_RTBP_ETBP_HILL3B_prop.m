% Script to compare propagated orbits from the three different  models:
% 1) RTBP
% 2) ETBP
% 3) Hill

% Constants and initial conditions (same as before)
earth_radius = 6378e3;
altitude = 300e3;
r0 = earth_radius + altitude;

G = 6.67430e-11;
M_earth = 5.972e24;
v0 = sqrt(G * M_earth / r0);

r0_nd = r0 / earth_radius;
v0_nd = v0 / (earth_radius / (24 * 3600));

% Set global parameters
global mu avg_e nu
mu = 0.012150582;
avg_e = 0.0549;
nu = 0;

% Initial states for each model
initial_state_rtbp = [mu - 1 + (mu/3)^(1/3); 0.01; 0.02; 0; -0.1; 0.02]; % Position and velocity % For RTBP (6 elements)
initial_state_etbp = [mu - 1 + (mu/3)^(1/3); 0.01; 0.02; 0; -0.1; 0.02; nu]; % For ETBP (7 elements = state vector + true anomaly)
initial_state_hill = [mu - 1 + (mu/3)^(1/3); 0.01; 0.02; 0; -0.1; 0.02]; % Position and velocity % For hill (6 elements)


% Time span for propagation (normalized to Earth-Moon system)
tspan = [0, 3]; % This range will likely need adjustment for full orbits

% Integrate the RTBP equations of motion
[t_rtbp, state_rtbp] = ode78(@HFEM_rtbp, tspan, initial_state_rtbp);

% Integrate the ETBP equations of motion
[t_etbp, state_etbp] = ode78(@HFEM_etbp, tspan, initial_state_etbp);

% Integrate Hill equations of motion
[t_hill, state_hill] = ode78(@HFEM_hill3B, tspan, initial_state_rtbp);

% Plot results
figure;
hold on;
plot3(state_rtbp(:,1), state_rtbp(:,2), state_rtbp(:,3), 'b', 'DisplayName', 'RTBP');
plot3(state_etbp(:,1), state_etbp(:,2), state_etbp(:,3), 'r', 'DisplayName', 'ETBP');
plot3(state_hill(:,1), state_hill(:,2), state_hill(:,3), 'g', 'DisplayName', 'Hill Model');

%plot3(mu, 0,0, 'DisplayName', 'Earth')
%plot3(mu-1, 0,0, 'DisplayName', 'Moon')
earth_moon_sphere('moon_dist', 'both') %fancy function to plot earth and moon...
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('RTBP Orbit', 'ETBP Orbit', 'Hill Orbit', 'Earth', 'Moon');
legend('show');
title('Comparison of RTBP, ETBP, and Hill Model Trajectories');
grid on;
hold off;


% Interpolate ETBP and Hill results to match the time points of RTBP results
state_etbp_interp = interp1(t_etbp, state_etbp(:,1:3), t_rtbp);
state_hill_interp = interp1(t_hill, state_hill(:,1:3), t_rtbp);
% Compute position differences for comparison
position_diff_rtbp_etbp = vecnorm(state_rtbp(:,1:3) - state_etbp_interp(:,1:3), 2, 2);
position_diff_rtbp_hill = vecnorm(state_rtbp(:,1:3) - state_hill_interp(:,1:3), 2, 2);

% Plot position differences over time
figure;
plot(t_rtbp, position_diff_rtbp_etbp, 'r', 'DisplayName', 'RTBP - ETBP');
hold on;
plot(t_rtbp, position_diff_rtbp_hill, 'g', 'DisplayName', 'RTBP - Hill');
xlabel('Time');
ylabel('Position Difference (nondimensional)');
legend('show');
title('Position Differences Between RTBP, ETBP, and Hill Model');
grid on;
hold off;

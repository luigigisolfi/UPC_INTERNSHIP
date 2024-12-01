function a = gravitational_acceleration_synodic(x)

global interpolators
global mu 

% -----------------------------------------------------------------%
% This function is used by full_force.m in order to retrieve the
% gravitational acceleration of a given body in the list 
% all_bodies = [PRIMARIES, BODIES] on the spacecraft
% -----------------------------------------------------------------%

BODY_str = regexprep(BODY, [{'\s+'}, {'-'}], '_');
% state_b = get_planet_state_vector(interpolators, BODY_str, t_inertial);

% r_b = state_b(1:3)

x1 = [mu,0,0];
x2 = [mu-1,0,0];

x1s = x(1:3) - x1(:);
x2s = x(1:3) - x2(:);

rho1s = (x1(1)-x(1))^2+(x1(2)-x(2))^2+(x1(3)-x(3))^2;
rho1s3 = rho1s*sqrt(rho1s);
rho2s = (x2(1)-x(1))^2+(x2(2)-x(2))^2+(x2(3)-x(3))^2;
rho2s3 = rho2s*sqrt(rho2s);

a1s = - (1-mu)*(x1s)/rho1s3;
a2s = - (mu)*(x2s)/rho2s3;

a = a1s +a2s;
% 
% fprintf('r_b-r_s')
% r_b(:)- r_s(:)
% r_sb2 = (r_s(1) - r_b(1))^2 + (r_s(2) - r_b(2))^2 + (r_s(3) - r_b(3))^2; %square of distance from spacecraft to body (in FRAME, OBSERVER)
% r_sb3 = r_sb2*sqrt(r_sb2);
% 
% a =  - GM_BODY*(x(:) - r_b(:))/r_sb3;


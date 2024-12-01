function a = gravitational_acceleration(t_inertial, x_inertial, BODY, GM_BODY)

global interpolators

% -----------------------------------------------------------------%
% This function is used by full_force.m in order to retrieve the
% gravitational acceleration of a given body in the list 
% all_bodies = [PRIMARIES, BODIES] on the spacecraft
% -----------------------------------------------------------------%

BODY_str = regexprep(BODY, [{'\s+'}, {'-'}], '_');
state_b = get_planet_state_vector(interpolators, BODY_str, t_inertial);

r_b = state_b(1:3);
r_s = x_inertial(1:3); %position of spacecraft

r_sb2 = (r_s(1) - r_b(1))^2 + (r_s(2) - r_b(2))^2 + (r_s(3) - r_b(3))^2; %square of distance from spacecraft to body (in FRAME, OBSERVER)
r_sb3 = r_sb2*sqrt(r_sb2);

a =  - GM_BODY*(r_s(:) - r_b(:))/r_sb3;


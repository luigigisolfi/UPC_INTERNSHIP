function xdot = HFEM_residual_etbp(t,x)

global PRIMARIES
global BODIES
global n_anomalistic
global mu
global eclipse_date_et
global avg_e
global nu

%here, t = rtbp_times. But interpolators are given in inertial physical
%coordinates. therefore, if we want to retrieve rs_rp, we need to write
%interpolators(inertial_t) where inertial_t = t*n_anomalistic = rtbp_times*n_anomalistic

inertial_t = (t + eclipse_date_et);

META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
inertial_state_body_1 = cspice_spkezr(PRIMARIES{1}, inertial_t, 'J2000', 'None', 'SUN');
inertial_state_body_1_plus = cspice_spkezr(PRIMARIES{1}, inertial_t+60, 'J2000', 'None', 'SUN');
inertial_state_body_1_plus_plus = cspice_spkezr(PRIMARIES{1}, inertial_t+120, 'J2000', 'None', 'SUN');
rp = inertial_state_body_1(1:3);
vp = inertial_state_body_1(4:6);
vp_plus = inertial_state_body_1_plus(4:6);
vp_plus_plus = inertial_state_body_1_plus_plus(4:6);
ap = (vp_plus - vp)/60;
ap_plus = (vp_plus_plus - vp_plus)/60;
oap = (ap_plus - ap)/60;

inertial_state_body_2 = cspice_spkezr(PRIMARIES{2}, inertial_t, 'J2000', 'None', 'SUN');
inertial_state_body_2_plus = cspice_spkezr(PRIMARIES{2}, inertial_t+60, 'J2000', 'None', 'SUN');
inertial_state_body_2_plus_plus = cspice_spkezr(PRIMARIES{2}, inertial_t+120, 'J2000', 'None', 'SUN');
rs = inertial_state_body_2(1:3);
vs = inertial_state_body_2(4:6);
vs_plus = inertial_state_body_2_plus(4:6);
vs_plus_plus = inertial_state_body_2_plus_plus(4:6);
as = (vs_plus - vs)/60;
as_plus = (vs_plus_plus - vs_plus)/60;
oas = (as_plus - as)/60;

%retrieving relative position, velocity and acceleration at time t
format long
rs_rp = rs - rp;
vs_vp = vs -vp;
as_ap = as - ap;
oas_oap = oas-oap;

SEb_pos = rp + mu*rs_rp;
SEb_vel = vp + mu*vs_vp;
SEb_acc = ap + mu*as_ap;

b = SEb_pos;
b_dot = SEb_vel;
b_ddot = SEb_acc;
k = norm(rs_rp);
k_dot = dot(rs_rp, vs_vp)/k;
k_ddot = (dot(vs_vp, vs_vp) + dot(rs_rp, as_ap))/k - (k_dot^2)/k;
h = norm(cross(rs_rp, vs_vp));
hp = dot(cross(rs_rp, as_ap), cross(rs_rp, vs_vp))/h;

%

C = construct_C(rs_rp, vs_vp);
C_dot = construct_C_dot(C, rs_rp, vs_vp, as_ap);
C_ddot = construct_C_ddot(C, C_dot, rs_rp, vs_vp, as_ap, oas_oap);

%
b1 = -b_ddot(1)/(n_anomalistic^2*k);
b2 = -b_ddot(2)/(n_anomalistic^2*k);
b3 = -b_ddot(3)/(n_anomalistic^2*k);
b4 = -2*k_dot/(n_anomalistic*k) + avg_e*sin(nu)/(2*sqrt(1+avg_e*cos(nu)));
b5 = 2*h/(n_anomalistic*k^2) - 2*(sqrt(1+avg_e*cos(nu)));
b6 = 2*k*as_ap(3)/(n_anomalistic*h);
b7 = -k_ddot/(n_anomalistic^2*k) + h^2/(n_anomalistic^2*k^4) -1;
b8 = - as_ap(3)/(n_anomalistic^2*k);
b9 = hp/(n_anomalistic^2*k^2);
b10 = -k_ddot/(n_anomalistic^2*k) + h^2/(n_anomalistic^2*k^4) + k^2*(as_ap(3))^2/(n_anomalistic^2*h^2)-1;
b11 = (3*h*k_dot - 2*k*hp)*(as_ap(3))/(n_anomalistic^2*h^2) + k*oas_oap(3)/(n_anomalistic^2*h);
b12 = - k_ddot/(n_anomalistic^2*k) + k^2*(as_ap(3))^2/(n_anomalistic^2*h^2) + avg_e*cos(nu);
b13 = 1;

first = [b1;b2;b3];
second = [b4,b5,0;-b5,b4,b6;0,-b6,b4];
third = [b7,b8,b9;-b8,b10,b11;b9,-b11,b12];

all_bodies = BODIES;
n_bodies = length(all_bodies);
xdot = zeros(6,1);

GM_1 = get_GM_body(PRIMARIES{1});
GM_2 = get_GM_body(PRIMARIES{2});

META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels

rtbp_pos_body_1 = (C\(b-rp)/k); 
rtbp_pos_body_2 = (C\(b-rs)/k);

x_s1 = x(1:3)-rtbp_pos_body_1;
x_s2 = x(1:3)-rtbp_pos_body_2;

rho_s13 = norm(x_s1)^3;
rho_s23 = norm(x_s2)^3;
synodic_acc_primaries = -(1-mu)*x_s1/rho_s13 - (mu)*x_s2/rho_s23; %primaries contribution to synodic acceleration

if isempty(BODIES) == 0

for i = 1:n_bodies % PARALLELIZED COMPUTATION ON ALL BODIES
    BODY = all_bodies{i};
    BODY_str = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    mu_body = get_GM_body(BODY)/(GM_1+GM_2);
    META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
    cspice_furnsh(META); %furnish kernel
    inertial_pos_body = cspice_spkpos(BODY_str, inertial_t, 'J2000', 'None', 'SUN');
    rtbp_pos_body =  (C\(b-inertial_pos_body)/k); 
    x_sb = x(1:3)-rtbp_pos_body;
    rho_sb3 = norm(x_sb)^3;
    synodic_acc_bodies = -mu_body*x_sb/rho_sb3; %primaries + bodies contribution to synodic acceleration
    acc_bodies(:, i) =synodic_acc_bodies;
    
end

Delta_Omega = sum(acc_bodies,2) + synodic_acc_primaries;

else
    Delta_Omega = synodic_acc_primaries;
    %fprintf('NO BODIES OTHER THAN THE PRIMARIES\n')
end
     %now compute the gravitational_acceleration_synodic(all bodies)
     %the function gravitational_acceleration gives acc in the inertial frame, so we need to divide the GMj contribution by (GM_1+GM_2) to get mu_j

xdot(1:3) = x(4:6);
xdot(4:6) = first + second*x(4:6) + third*x(1:3) + b13*Delta_Omega;
xdot(7) = sqrt(1+avg_e*cos(nu));

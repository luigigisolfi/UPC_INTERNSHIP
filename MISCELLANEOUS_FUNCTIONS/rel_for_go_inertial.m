function rel_for_go_inertial(t)
global PRIMARIES
global BODIES
global interpolators
global n_anomalistic
global mu
global eclipse_date_et

%----------------------------------------------------------------------------
% SORRY! I AM NOT SURE WHAT THIS FUNCTION WAS FOR. ANYWAY, IT IS PROBABLY
% NOT NEEDED, AND HAS BEEN SUPERSEEDED BY SOME COOLER ONES...
%-----------------------------------------------------------------------------

%here, t = rtbp_times. But interpolators are given in inertial physical
%coordinates. therefore, if we want to retrieve rs_rp, we need to write
%interpolators(inertial_t) where inertial_t = t*n_anomalistic = rtbp_times*n_anomalistic

inertial_t = (t+ eclipse_date_et);
% Evaluate the splines for the primaries
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
    interpolated_p = ppval(interpolators.(primary_str).spline{dim}, inertial_t);
    interpolated_s = ppval(interpolators.(secondary_str).spline{dim}, inertial_t);
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

%retrieving relative acceleration at time t
rs_rp = rp-rs;
vs_vp = vp-vs;
as_ap = ap-as;
oas_oap = oap-oas;

SEb_pos = rp - mu*rs_rp;
SEb_vel = vp - mu*vs_vp;
SEb_acc = ap - mu*as_ap;

b = SEb_pos;
b_dot = SEb_vel;
b_ddot = SEb_acc;
k = norm(rs_rp);
k_dot = dot(rs_rp, vs_vp)/k;
k_ddot = (dot(vs_vp, vs_vp) + dot(rs_rp, as_ap))/k - (k_dot^2)/k;
h = norm(cross(rs_rp, vs_vp));
hp = dot(cross(rs_rp, as_ap), cross(rs_rp, vs_vp))/h;
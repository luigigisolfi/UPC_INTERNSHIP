function C_ddot = construct_C_ddot(C, C_dot, rs_rp, vs_vp, as_ap, oas_oap)

%-------------------------------------------------------------------------%
% This function builds the C'' matrix For acceleration conversions 
% [see Luigi Gisolfi's internship @UPC report: 
% https://www.overleaf.com/project/6656db075c745bb769c47f7c 
% for full derivation]

C1 = C(:,1);
C2 = C(:,2);
C3 = C(:,3);

C1_dot = C_dot(:,1);
C2_dot = C_dot(:,2);
C3_dot = C_dot(:,3);

h = norm(cross(rs_rp, vs_vp));
hp = dot(cross(rs_rp, as_ap), cross(rs_rp, vs_vp))/h;

k = norm(rs_rp);
k_dot = dot(rs_rp, vs_vp)/k;
k_ddot = (dot(vs_vp, vs_vp) + dot(rs_rp, as_ap))/k - (k_dot^2)/k;

C1_ddot = as_ap/k - 2*k_dot*vs_vp/k^2 + (2*k_dot^2 - k*k_ddot)*rs_rp/k^3;
etemp = cross(rs_rp, as_ap) + cross(rs_rp, oas_oap);
C3_ddot = etemp/h - 2*hp*(cross(rs_rp, as_ap))/h^2;
C2_ddot = cross(C3_ddot, C1) + 2*cross(C3_dot, C1) + cross(C3, C1_ddot);
C_ddot = [C1_ddot, C2_ddot, C3_ddot];

end
%-------------------------------------------------------------------------%

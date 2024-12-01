function C_dot = construct_C_dot(C, rs_rp, vs_vp, as_ap)

%-------------------------------------------------------------------------%
% This function builds the C' matrix For acceleration conversions 
% [see Luigi Gisolfi's internship @UPC report: 
% https://www.overleaf.com/project/6656db075c745bb769c47f7c 
% for full derivation]

C1 = C(:,1);
C2 = C(:,2);
C3 = C(:,3);


h = norm(cross(rs_rp, vs_vp));
hp = dot(cross(rs_rp, as_ap), cross(rs_rp, vs_vp))/h;

k = norm(rs_rp);
k_dot = dot(rs_rp, vs_vp)/k;

C1_dot = (vs_vp - k_dot*C1)/k;
C3_dot = (cross(rs_rp, as_ap) - hp*C3)/h;
C2_dot = cross(C3_dot, C1) + cross(C3, C1_dot);
C_dot = [C1_dot, C2_dot, C3_dot];

end
%-------------------------------------------------------------------------%

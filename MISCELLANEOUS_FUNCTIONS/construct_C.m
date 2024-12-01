function C = construct_C(rs_rp, vs_vp)

%-------------------------------------------------------------------------%
% This function builds the C matrix For acceleration conversions 
% [see Luigi Gisolfi's internship @UPC report: 
% https://www.overleaf.com/project/6656db075c745bb769c47f7c 
% for full derivation]

k = norm(rs_rp);
C1 = rs_rp/k;
C3 = cross(rs_rp, vs_vp)/norm(cross(rs_rp, vs_vp));
C2 = cross(C3, C1);
C2 = C2;
C = [C1, C2, C3];
end
%-------------------------------------------------------------------------%

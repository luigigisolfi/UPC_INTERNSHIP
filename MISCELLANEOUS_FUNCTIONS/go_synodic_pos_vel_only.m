function [rtbp_pos_spacecraft, rtbp_vel_spacecraft] = go_synodic_pos_vel_only(rel_pos, rel_vel, rel_acc, SEb_pos, SEb_vel, inertial_state, n)

% This function transforms synodic quantities into inertial ones, as done
% in "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points"
% Inputs: 
% rel_pos, rel_vel, rel_acc, over_acc (relative, of the two primaries)
% SEb_pos, SEb_vel, SEb_acc (barycenter)
% inertial_pos, inertial_vel, inertial_acc (rtbp coordinates to be converted)
% n (mean anomaly value)
%
% Output:
% rtbp_pos_spacecraft, rtbp_vel_spacecraft

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nPerforming conversion from synodic (adimensional) system to inertial (physical) system...\n')
rtbp_pos_spacecraft = zeros(3,length(rel_pos));
rtbp_vel_spacecraft = zeros(3,length(rel_pos));

for i = 1:length(rel_pos)
    rs_rp = rel_pos(:, i); 
    k = norm(rs_rp);
    vs_vp = rel_vel(:, i); 
    ap_as = rel_acc(:,i);
    

    C = construct_C(rs_rp, vs_vp);
    C_dot = construct_C_dot(C, rs_rp, vs_vp, ap_as);

    b = SEb_pos(1:3, i);
    b_dot = SEb_vel(:,i);
    k_dot = dot(rs_rp, vs_vp)/k;

 
    %Apply the formula to conversion between two reference systems 
    %(see "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points", pp. 137-138)
    rtbp_pos_spacecraft(:,i) = C\(inertial_state(i,1:3).'-b)/k;  
    rtbp_vel_spacecraft(:,i) = (C*k*n)\(inertial_state(i, 4:6).' - b_dot - k*C_dot*rtbp_pos_spacecraft(:,i) - k_dot*C*rtbp_pos_spacecraft(:,i)); 
end
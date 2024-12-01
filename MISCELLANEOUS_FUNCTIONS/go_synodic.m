function [rtbp_pos_spacecraft, rtbp_vel_spacecraft, rtbp_acc_spacecraft] = go_synodic(rel_pos, rel_vel, rel_acc, over_acc, SEb_pos, SEb_vel, SEb_acc, inertial_pos, inertial_vel, inertial_acc, n)

% This function transforms synodic quantities into inertial ones, as done
% in "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points"
% Inputs: 
% rel_pos, rel_vel, rel_acc, over_acc (relative, of the two primaries)
% SEb_pos, SEb_vel, SEb_acc (barycenter)
% inertial_pos, inertial_vel, inertial_acc (rtbp coordinates to be converted)
% n (mean anomaly value)
%
% Output:
% rtbp_pos_spacecraft, rtbp_vel_spacecraft, rtbp_acc_spacecraft

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nPerforming conversion from synodic (adimensional) system to inertial (physical) system...\n')
rtbp_pos_spacecraft = zeros(3,length(rel_pos));
rtbp_vel_spacecraft = zeros(3,length(rel_pos));
rtbp_acc_spacecraft = zeros(3, length(rel_pos));

check_inertial_pos = zeros(3, length(rel_pos));
check_inertial_vel = zeros(3, length(rel_pos));

for i = 1:length(rel_pos)
    rs_rp = rel_pos(:, i); 
    k = norm(rs_rp);
    vs_vp = rel_vel(:, i); 
    as_ap = rel_acc(:,i);
    
    oas_oap = over_acc(:,i);


    C = construct_C(rs_rp, vs_vp);
    C_dot = construct_C_dot(C, rs_rp, vs_vp, as_ap);
    C_ddot = construct_C_ddot(C, C_dot, rs_rp, vs_vp, as_ap, oas_oap);

    b = SEb_pos(1:3, i);
    b_dot = SEb_vel(:,i);
    b_ddot = SEb_acc(:,i);
    k_dot = dot(rs_rp, vs_vp)/k;
    k_ddot = (dot(vs_vp, vs_vp) + dot(rs_rp, as_ap) - k_dot^2)/k;
    

    %Apply the formula to conversion between two reference systems 
    %(see "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points", pp. 137-138)
    rtbp_pos_spacecraft(:,i) =  C\(inertial_pos(:, i)-b)/k; 

    rtbp_vel_spacecraft(:,i) = (C*k*n)\(inertial_vel(:,i) - b_dot - k*C_dot*rtbp_pos_spacecraft(:,i) - k_dot*C*rtbp_pos_spacecraft(:,i)); 

    a = rtbp_pos_spacecraft(:,i);
    a_dot = rtbp_vel_spacecraft(:,i);
    R_ddot = inertial_acc(:,i);
    rtbp_acc_spacecraft(:,i) =(C*k*n^2)\(R_ddot - b_ddot - (k_ddot*C + 2*k_dot*C_dot+ k*C_ddot)*a - (2*k_dot*C + 2*k*C_dot)*a_dot*n);

    %--------------------------------------------------------------------------------------------%   
    % These two arrays are needed for the check of positions and velocities
    % converted back to rtbp coordinates (if requested)
    %check_inertial_pos(:,i) = k*C*rtbp_pos_spacecraft(:,i)+ b ;
    %check_inertial_vel(:,i) =b_dot + k_dot*C*rtbp_pos_spacecraft(:,i) + k*(C_dot*rtbp_pos_spacecraft(:,i)+ n*C*rtbp_vel_spacecraft(:,i));
   
end

%fprintf('Conversion: Done.\n')
%--------------------------------------------------------------------------------------------%
% Perform a quick check that the obtained inertial positions and velocities
% can be convertedback into the original rtbp positions and velocities and
% they are the same if converted back to rtbp coordinates

%fprintf('-----------------------------------------------------------\n')
%fprintf('Function: go_inertial\nChecking whether the conversion has been successfull...\n')
%array_to_check = check_inertial_pos - inertial_pos;
%check = 0;
%for i = 1:3
%    for j = 1:length(array_to_check)
%       if array_to_check(i,j) > 1e-10
%           check = 1 ;
%       end
%    end
%end

%if check == 1
%    error('Attention. Conversion check has not been successfull. Aborting...\n')
%    fprintf('-----------------------------------------------------------\n')
%else
%    fprintf('All good.\n')
%end
%--------------------------------------------------------------------------------------------%

function [inertial_pos_spacecraft, inertial_vel_spacecraft, inertial_acc_spacecraft] = go_inertial(rel_pos, rel_vel, rel_acc, over_acc, SEb_pos, SEb_vel, SEb_acc, rtbp_pos, rtbp_vel, rtbp_acc, n)
global mu
global inertial_pos_earth
global inertial_pos_moon


% This function transforms synodic quantities into inertial ones, as done
% in "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points"
% Inputs: 
% rel_pos, rel_vel, rel_acc, over_acc (relative, of the two primaries)
% SEb_pos, SEb_vel, SEb_acc (barycenter)
% rtbp_pos, rtbp_vel, rtbp_acc (rtbp coordinates to be converted)
% n (mean anomaly value)
%
% Output:
% inertial_pos_spacecraft, inertial_vel_spacecraft, inertial_acc_spacecraft

%NOTE:
% The function also checks that the conversion has been successfull by
% re-converting back to rtbp coordinates and making sure the difference is
% within a given threshold

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nPerforming conversion from synodic (adimensional) system to inertial (physical) system...\n')
inertial_pos_spacecraft = zeros(3,length(rel_pos));
inertial_vel_spacecraft = zeros(3,length(rel_pos));
inertial_acc_spacecraft = zeros(3, length(rel_pos));
inertial_pos_earth =zeros(3, length(rel_pos));
inertial_pos_moon =zeros(3, length(rel_pos));
check_rtbp_pos = zeros(3, length(rel_pos));
check_rtbp_vel = zeros(3, length(rel_pos));

for i = 1:length(rel_pos)
    rs_rp = rel_pos(:,i);
    k = vecnorm(rs_rp);
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
    inertial_pos_spacecraft(:,i) = k*C*rtbp_pos(:,i)+ b ;
    inertial_vel_spacecraft(:,i) = b_dot + k_dot*C*rtbp_pos(:,i) + k*(C_dot*rtbp_pos(:,i)+ n*C*rtbp_vel(:,i));
    inertial_pos_earth(:,i) = k*C*[mu,0,0].'+ b;
    inertial_pos_moon(:,i) = k*C*[mu-1,0,0].'+ b;
    a = rtbp_pos(:,i);
    a_dot = rtbp_vel(:,i);
    a_ddot = rtbp_acc(:,i);
    inertial_acc_spacecraft(:,i) = b_ddot + (k_ddot*C + 2*k_dot*C_dot+ k*C_ddot)*a + (2*k_dot*C + 2*k*C_dot)*a_dot*n + k*C*n^2*a_ddot;

    %--------------------------------------------------------------------------------------------%   
    % These two arrays are needed for the check of positions and velocities
    % converted back to rtbp coordinates (if requested)
    check_rtbp_pos(:,i) = C\(inertial_pos_spacecraft(:,i)-b)/k; 
    check_rtbp_vel(:,i) =(C*k*n)\(inertial_vel_spacecraft(:,i) - b_dot - k*C_dot*rtbp_pos(:,i) - k_dot*C*rtbp_pos(:,i)); %store velocities in structure
   
end

fprintf('Conversion: Done.\n')
%--------------------------------------------------------------------------------------------%
% Perform a quick check that the obtained inertial positions and velocities
% can be convertedback into the original rtbp positions and velocities and
% they are the same if converted back to rtbp coordinates

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nChecking whether the conversion has been successfull...\n')
array_to_check = check_rtbp_pos - rtbp_pos;
check = 0;
for i = 1:3
    for j = 1:length(array_to_check)
       if array_to_check(i,j) > 1e-10
           check = 1 ;
       end
    end
end

if check == 1
    error('Attention. Conversion check has not been successfull. Aborting...\n')
    fprintf('-----------------------------------------------------------\n')
else
    fprintf('All good.\n')
end
%--------------------------------------------------------------------------------------------%

function [rtbp_pos_spacecraft] = go_synodic_pos_only(rel_pos, rel_vel, SEb_pos, inertial_state, n)

% This function transforms synodic quantities into inertial ones, as done
% in "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points"
% Inputs: 
% rel_pos, rel_vel, rel_acc, over_acc (relative, of the two primaries)
% SEb_pos, SEb_vel, SEb_acc (barycenter)
% inertial_pos, inertial_vel, inertial_acc (rtbp coordinates to be converted)
% n (mean anomaly value)
%
% Output:
% rtbp_pos_spacecraft

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nPerforming conversion from synodic (adimensional) system to inertial (physical) system...\n')
rtbp_pos_spacecraft = zeros(3,length(rel_pos));

check_inertial_pos = zeros(3, length(rel_pos));
check_inertial_vel = zeros(3, length(rel_pos));

for i = 1:length(rel_pos)
    rs_rp = rel_pos(:, i); 
    k = norm(rs_rp);
    vs_vp = rel_vel(:, i); 

    C = construct_C(rs_rp, vs_vp);

    b = SEb_pos(1:3, i);

    
    %Apply the formula to conversion between two reference systems 
    %(see "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points", pp. 137-138)

    rtbp_pos_spacecraft(:,i) =  C\(inertial_state(i,1:3).'-b)/k; 

    %--------------------------------------------------------------------------------------------%   
    % These two arrays are needed for the check of positions and velocities
    % converted back to rtbp coordinates (if requested)
    check_inertial_pos(:,i) = k*C*rtbp_pos_spacecraft(:,i)+ b ;
   
end

%fprintf('Conversion: Done.\n')
%--------------------------------------------------------------------------------------------%
% Perform a quick check that the obtained inertial positions and velocities
% can be convertedback into the original rtbp positions and velocities and
% they are the same if converted back to rtbp coordinates

fprintf('-----------------------------------------------------------\n')
fprintf('Function: go_inertial\nChecking whether the conversion has been successfull...\n')
array_to_check = check_inertial_pos - inertial_state(:,1:3).';
check = 0;
for i = 1:3
   for j = 1:length(array_to_check)
      if array_to_check(i,j) > 1e-7
          array_to_check(i,j)
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

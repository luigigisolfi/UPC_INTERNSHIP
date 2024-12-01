function [inertial_state_primaries, inertial_state_bodies, interpolators] = get_ephemeris(times, PRIMARIES, BODIES, FRAME, OBSERVER)
fprintf('Function: get_ephemeris\nRetrieving INERTIAL state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE)...\n')
inertial_state_primaries = struct();

for j = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{j};
    [inertial_state_primary, ~] = cspice_spkezr(PRIMARY, times, FRAME, 'NONE', OBSERVER); %State BODY wrt EM Barycenter
    % Compute the Relative Accelerations (needed for velocity conversion)
    acc = zeros(3,length(inertial_state_primary));
    acc(1,1:end-1) = diff(inertial_state_primary(4,:))./diff(times);
    acc(2,1:end-1) = diff(inertial_state_primary(5,:))./diff(times);
    acc(3,1:end-1) = diff(inertial_state_primary(6,:))./diff(times);
    acc(:,1) = acc(:,2);
    acc(:,end) = acc(:,end-1);
    % Compute the Relative Over Accelerations (needed for acceleration conversion)
    over_acc = zeros(3,length(inertial_state_primary));
    over_acc(1,1:end-1) = diff(acc(1,:))./diff(times);
    over_acc(2,1:end-1) = diff(acc(2,:))./diff(times);
    over_acc(3,1:end-1) = diff(acc(3,:))./diff(times);
    over_acc(:,1) = over_acc(:,2);
    over_acc(:,end) = over_acc(:,end-1);
    
    % Check if the string contains a hyphen% Replace all backspaces and hyphens with underscores
    PRIMARY_str{j} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    inertial_state_primaries.(PRIMARY_str{j}).position = inertial_state_primary(1:3, :);
    inertial_state_primaries.(PRIMARY_str{j}).velocity = inertial_state_primary(4:6, :);
    inertial_state_primaries.(PRIMARY_str{j}).acceleration = acc(1:3, :);
    inertial_state_primaries.(PRIMARY_str{j}).overacceleration = over_acc(1:3, :);
    
end

fprintf('Retrieved.\n\n')
fprintf('Retrieving INERTIAL state vectors of the BODIES in the ECLIPJ2000 frame (from SPICE)...\n')
inertial_state_bodies = struct();
for j = 1:length(BODIES)
    BODY = BODIES{j};
    [inertial_state_body, ~] = cspice_spkezr(BODY, times, FRAME, 'NONE', OBSERVER); %State BODY wrt SS Barycenter
    % Compute the Relative Accelerations (needed for velocity conversion)
    acc_b = zeros(3,length(inertial_state_body));
    acc_b(1,1:end-1) = diff(inertial_state_body(4,:))./diff(times);
    acc_b(2,1:end-1) = diff(inertial_state_body(5,:))./diff(times);
    acc_b(3,1:end-1) = diff(inertial_state_body(6,:))./diff(times);
    acc_b(:,1) = acc_b(:,2);
    acc_b(:,end) = acc_b(:,end-1);

    over_acc_b = zeros(3,length(inertial_state_body));
    over_acc_b(1,1:end-1) = diff(acc_b(1,:))./diff(times);
    over_acc_b(2,1:end-1) = diff(acc_b(2,:))./diff(times);
    over_acc_b(3,1:end-1) = diff(acc_b(3,:))./diff(times);
    over_acc_b(:,1) = over_acc_b(:,2);
    over_acc_b(:,end) = over_acc_b(:,end-1);

    % Check if the string contains a hyphen% Replace all backspaces and hyphens with underscores
    BODY_str{j} = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    inertial_state_bodies.(BODY_str{j}).position = inertial_state_body(1:3, :);
    inertial_state_bodies.(BODY_str{j}).velocity = inertial_state_body(4:6, :);
    inertial_state_bodies.(BODY_str{j}).acceleration = acc_b(1:3, :);
    inertial_state_bodies.(BODY_str{j}).overacceleration = over_acc_b(1:3, :);
end

fprintf('Retrieved.\n\n')
%-------- Create spline interpolators for state and vector field of PRIMARIES ---------- %

fprintf('Creating cubic Chebychev interpolators for PRIMARIES and BODIES...\n')
% Store merged position and velocity data for each primary
for j = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{j};
    PRIMARY_str{j} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    inertial_state_primaries_array.(PRIMARY_str{j}) = [inertial_state_primaries.(PRIMARY_str{j}).position; inertial_state_primaries.(PRIMARY_str{j}).velocity; inertial_state_primaries.(PRIMARY_str{j}).acceleration; inertial_state_primaries.(PRIMARY_str{j}).overacceleration];
end

interpolator_primaries = struct();
for i = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{i};
    PRIMARY_str{i} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    interpolator_primaries.(PRIMARY_str{i} ).x = times;
    interpolator_primaries.(PRIMARY_str{i} ).y = inertial_state_primaries_array.(PRIMARY_str{i});
    interpolator_primaries.(PRIMARY_str{i} ).spline = cell(1,12);  % For 6 dimensions (3 for position, 3 for velocity)
    for dim = 1:12
        interpolator_primaries.(PRIMARY_str{i} ).spline{dim} = spline(times, inertial_state_primaries_array.(PRIMARY_str{i} )(dim,:));
    end
end
%-------- Create spline interpolators for state and vector field of BODIES ---------- %
 % Store merged position and velocity data for each primary
for j = 1:length(BODIES)
    BODY = BODIES{j};
 
    PRIMARY_str{j} = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    inertial_state_bodies_array.(BODY_str{j}) = [inertial_state_bodies.(BODY_str{j}).position; inertial_state_bodies.(BODY_str{j}).velocity; inertial_state_bodies.(BODY_str{j}).acceleration; inertial_state_bodies.(BODY_str{j}).overacceleration];
end

for i = 1:length(BODIES)
    BODY = BODIES{i};
    BODY_str{i} = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    interpolator_bodies.(BODY_str{i} ).x = times;
    interpolator_bodies.(BODY_str{i} ).y = inertial_state_bodies_array.(BODY_str{i});
    interpolator_bodies.(BODY_str{i} ).spline = cell(1,12);  % For 6 dimensions (3 for position, 3 for velocity)
    for dim = 1:12
        interpolator_bodies.(BODY_str{i} ).spline{dim} = spline(times, inertial_state_bodies_array.(BODY_str{i} )(dim,:));
    end
end

%------ Combine interpolators for primaries and bodies into one struct ------ %
interpolators = struct();

% Add interpolators for primaries
for i = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{i};
    PRIMARY_str{i} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    interpolators.(PRIMARY_str{i}) = interpolator_primaries.(PRIMARY_str{i});
end

% Add interpolators for bodies
for i = 1:length(BODIES)
    BODY_str{i} = regexprep(BODY_str{i}, [{'\s+'}, {'-'}], '_');
    interpolators.(BODY_str{i}) = interpolator_bodies.(BODY_str{i});
end
fprintf('Created.\n')

end
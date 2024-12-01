function [rtbp_times, rtbp_state_primaries, body_1, body_2, object, rtbp_state_object, interpolators_object] = convert_object(ti, tf, FRAME, OBSERVER, PRIMARIES, L, object)

% Function: spice_to_rtbp 
% This function:
% 1) Retrieves SPICE ephemeris data for a given body, in the FRAME "ECLIPJ2000", wtr to an OBSERVER (chosen by the user)
% 2) Performs a change of coordinates from the SPICE frame coordinates (and units) to the RTBP coordinates (and units)

% Inputs are:
% ti, tf: user chosen initial and final epochs (in TDB format)
% FRAME: Initial Reference Frame
% OBSERVER: Origin of the Frame
% PRIMARIES: primaries to be Considered
% object: Additional Solar System body/object for which we want to retrieve ephemeris (i.e. Venus, Mercury, Mars, ecc...)

% Outputs are:
% rtbp_times: rtbp-like adimensional de40X epochs in normalized units; 
% rtbp_state_primaries: Primaries' rtbp-like adimensional state vectors in normalized units;
% body_1, body_2: Name of PRIMARIES considered;
% object: name of object considered
% rtbp_state_object: structure with fields corresponding to each BODY. Each field is a matrix of positions and velocities at each time;
% interpolators: the interpolated SPICE data, to be called in the full_force.m function to retrieve the position of object at each desired propagation epoch;

% ----- DEFINE GLOBAL VARIABLES ------- %
global mu
%---------------------------------------

% ----- Load required Kernels and ephemeris file ------ %
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
filename = '/Users/luigigisolfi/Documents/mcodes/mcodes/kernels/spk/de405.bsp';
% ------------------------------------------------------ %

% --------------- Retrieve object ids and corresponding names ------------- %
body_ids = cspice_spkobj(filename, 15); %Find the set of ID codes of all objects in a specified SPK file
wanted_ids = body_ids(1:9);
body_names = cell(size(wanted_ids));
% Loop over each body ID and get the corresponding name
for i = 1:length(body_ids)
    try
        % Try to get the body name using the ID
        body_names{i} = cspice_bodc2n(body_ids(i));
    catch ME
        % Handle the case where the body ID does not have a corresponding name
        if strcmp(ME.identifier, 'SPICE(NOFRAME)')
            body_names{i} = ['Unknown ID: ', num2str(body_ids(i))];
        else
            rethrow(ME); % Re-throw the error if it's not the expected one
        end
    end
end
% ------------------------------------------------------------ %

% -------------------------- Define the start and stop of the simulation --------------------- %
% Initialize the coverage window
cover = cspice_spkcov(filename, body_ids(1), 1);

% Extract the start time from the coverage window
start_time_et = cover(1);% The first element is the start time in ephemeris time (ET)
end_time_et = cover(end);

wanted_start_time_et = cspice_str2et(ti);
wanted_end_time_et = cspice_str2et(tf);

% Convert the start time from ET to a human-readable calendar format
start_time_str = cspice_et2utc(start_time_et, 'C', 0);
end_time_str = cspice_et2utc(end_time_et, 'C', 0);

fprintf('SPICE Ephemeris Start Epoch: %t\n', start_time_str)
fprintf('SPICE Ephemeris End Epoch: %t\n', end_time_str)

wanted_start_time_str = cspice_et2utc(wanted_start_time_et, 'C', 0);
wanted_end_time_str = cspice_et2utc(wanted_end_time_et, 'C', 0);

if wanted_start_time_et > end_time_et
     error('Start Simulation Epoch: %s is bigger than end of de405 coverage.\n Aborting...\n', wanted_start_time_str);
elseif wanted_end_time_et < start_time_et
     error('End Simulation Epoch: %s is smaller than start of de405 coverage.\n Aborting...\n', wanted_start_time_str);

elseif wanted_start_time_et > start_time_et && wanted_end_time_et < end_time_et
    fprintf('Start Simulation Epoch: %s\n', wanted_start_time_str);
    fprintf('End Simulation Epoch: %s\n', wanted_end_time_str);
elseif wanted_start_time_et > start_time_et && wanted_end_time_et > end_time_et
    error('End Simulation Epoch: %s is too big\n Aborting...', wanted_end_time_str);
elseif wanted_start_time_et < start_time_et && wanted_end_time_et < end_time_et
    error('Start Simulation Epoch: %s is too small\n Aborting...', wanted_start_time_str);
else
    error('Both start and end Simulation epochs fall outside of the de405 time coverage\n Aborting...\n');
end

STEPS = 10000;
times = linspace(wanted_start_time_et, wanted_end_time_et, STEPS);
% -----------------------------------------------------------------

%-------- Conversion from dimensional to adimensional time epochs ---------------
n = get_n(PRIMARIES, L); %L is the distance between the PRIMARIES in the assumption of circular orbit
rtbp_times = (times(:))*n;
%-------- Conversion from dimensional to adimensional time epochs ---------------

% ------------------------ RTBP variables ---------------------------%
[mu, body_1, body_2] = get_mu(PRIMARIES);
GM_1 = get_GM_body(body_1);
GM_2 = get_GM_body(body_2);
%---------------------------------------------------------------------%


%--------- Retrieve state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE) ---------%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
state_primaries = struct();
for j = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{j};
    [state_primary, ~] = cspice_spkezr(PRIMARY, times, FRAME, 'NONE', OBSERVER); %State BODY wrt SS Barycenter
    % Check if the string contains a hyphen% Replace all backspaces and hyphens with underscores
    PRIMARY_str{j} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    state_primaries.(PRIMARY_str{j}).position = state_primary(1:3, :);
    state_primaries.(PRIMARY_str{j}).velocity = state_primary(4:6, :);
end

SEb_pos = (GM_1 * state_primaries.(PRIMARY_str{1}).position + GM_2 * state_primaries.(PRIMARY_str{2}).position) / (GM_2 + GM_1);
SEb_vel = (GM_1 * state_primaries.(PRIMARY_str{1}).velocity + GM_2 * state_primaries.(PRIMARY_str{2}).velocity) / (GM_2 + GM_1);
%-----------------------------------------------------------------------------------------------------%

%--------- Retrieve state vectors of considered object in the ECLIPJ2000 frame (from SPICE) ---------%
state_object = struct();
for j = 1:length(object)
    BODY = object{j};
    BODY_str{j} = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    [state_body, ~] = cspice_spkezr(BODY, times, FRAME, 'NONE', OBSERVER); %State BODY wrt SS Barycenter
    state_object.(BODY_str{j} ).position = state_body(1:3, :);
    state_object.(BODY_str{j} ).velocity = state_body(4:6, :);
end
%-----------------------------------------------------------------------------------------------------%

% --------------------- GO FROM ECLIPTIC TO SYNODIC REFERENCE SYSTEM  -------------------------
pos_1 = state_primaries.(PRIMARY_str{1}).position;
pos_2 = state_primaries.(PRIMARY_str{2}).position;
rel_pos = pos_1 - pos_2;
vel_1 = state_primaries.(PRIMARY_str{1}).velocity;
vel_2 = state_primaries.(PRIMARY_str{2}).velocity;
rel_vel = vel_1 - vel_2;
rel_acc = zeros(3,STEPS);


rel_acc(1,1:end-1) = diff(rel_vel(1,:))./diff(times);
rel_acc(2,1:end-1) = diff(rel_vel(1,:))./diff(times);
rel_acc(3,1:end-1) = diff(rel_vel(1,:))./diff(times);
rel_acc(:,1) = rel_acc(:,2);
rel_acc(:,end) = rel_acc(:,end-1);

rtbp_pos_primaries = struct();
rtbp_vel_primaries = struct();

rtbp_pos_object = struct();
rtbp_vel_object = struct();

%Create the orthogonal matrix C
for i = 1:length(rel_pos) 
    rp_rs = rel_pos(:, i); 
    rp_rs_norm = norm(rp_rs);
    vp_vs = rel_vel(:, i); 
    ap_as = rel_acc(:,i);
    C1 = rp_rs/rp_rs_norm;
    C3 = cross(rp_rs, vp_vs)/norm(cross(rp_rs, vp_vs));
    C2 = cross(C3, C1);
    C = [C1, C2, C3];
    b = SEb_pos(1:3, i);
    b_dot = SEb_vel(:,i);
    k_dot = dot(rp_rs, vp_vs)/rp_rs_norm;
    
    C_dot_1 = (rp_rs_norm*(vp_vs) - k_dot*(rp_rs))/(rp_rs_norm^2);
    A = cross(rp_rs, ap_as)/norm(cross(rp_rs, vp_vs));
    B = dot(C3, (cross(rp_rs, ap_as)))/norm(cross(rp_rs,vp_vs)); %Correct version of eq 1.23 in Treball
    C_dot_3 = A - B;
    C_dot_2 = cross(C_dot_3, C1) + cross(C3, C_dot_1);
    C_dot = [C_dot_1, C_dot_2, C_dot_3];

    %Apply the formula to conversion between two reference systems 
    %(see "Jorba, Simo, Masdemont, Gomez, Dynamics and Mission Design Near Libration Points", pp. 137-138)

    for p = 1:length(PRIMARIES)
        PRIMARY = PRIMARIES{p};
        PRIMARY_str{j} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
        state_primary = state_primaries.(PRIMARY_str{j});

        rtbp_pos_primaries.(PRIMARY_str{j})(:,i) = C\(state_primary.position(:, i)-b)/rp_rs_norm; %store positions in structure
        
        (state_primary.velocity(:,i) - b_dot - rp_rs_norm*C_dot*rtbp_pos_primaries.(PRIMARY_str{j})(:,i) - k_dot*C*rtbp_pos_primaries.(PRIMARY_str{j})(:,i)); %store velocities in structure
        
        rtbp_vel_primaries.(PRIMARY_str{j})(:,i) = (C*rp_rs_norm*n)\(state_primary.velocity(:,i) - b_dot - rp_rs_norm*C_dot*rtbp_pos_primaries.(PRIMARY_str{j})(:,i) - k_dot*C*rtbp_pos_primaries.(PRIMARY_str{j})(:,i)); %store velocities in structure
    end
    
    for j = 1:length(object)
        BODY = object{j};
        BODY_str{j} = regexprep(BODY, [{'\s+'}, {'-'}], '_');
        state_body = state_object.(BODY_str{j});
        rtbp_pos_object.(BODY_str{j})(:,i) = C\(state_body.position(:, i)-b)/rp_rs_norm; %store positions in structure
        rtbp_vel_object.(BODY_str{j})(:,i) = (C*rp_rs_norm*n)\(state_body.velocity(:,i) - b_dot - rp_rs_norm*C_dot*rtbp_pos_object.(BODY_str{j})(:,i) - k_dot*C*rtbp_pos_object.(BODY_str{j})(:,i)); %store velocities in structure
    end
end


%------------- Print an output providing some context about the system ---------- %
fprintf('Initial Epoch (Normalized Units): %e\n', rtbp_times(1))
fprintf('Final Epoch (Normalized Units): %e\n', rtbp_times(end))
%--------------------------------------------------------------------------------%

% Store merged position and velocity data for each body
rtbp_state_object = struct();
for j = 1:length(object)
    BODY = object{j};
    BODY_str{j} = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    rtbp_state_object.(BODY_str{j}) = [rtbp_pos_object.(BODY_str{j}); rtbp_vel_object.(BODY_str{j})];
end

% Create spline interpolators for each planet's state
interpolators_object = struct();
for i = 1:length(object)
    BODY = object{i};
    BODY_str{j} = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    interpolators_object.(BODY_str{j}).x = rtbp_times;
    interpolators_object.(BODY_str{j}).y = rtbp_state_object.(BODY_str{j});
    interpolators_object.(BODY_str{j}).spline = cell(1,6);  % For 6 dimensions (3 for position, 3 for velocity)
    for dim = 1:6
        interpolators_object.(BODY_str{j}).spline{dim} = spline(rtbp_times, rtbp_state_object.(BODY_str{j})(dim,:));
    end
end

% Store merged position and velocity data for each body
rtbp_state_primaries= struct();
for j = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{j};
 
    PRIMARY_str{j} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    rtbp_state_primaries.(PRIMARY_str{j}) = [rtbp_pos_primaries.(PRIMARY_str{j}); rtbp_vel_primaries.(PRIMARY_str{j})];
end

% Create spline interpolators for each planet's state
interpolators_primaries = struct();
for i = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{i};
    PRIMARY_str{i} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    interpolators_primaries.(PRIMARY_str{i} ).x = rtbp_times;
    interpolators_primaries.(PRIMARY_str{i} ).y = rtbp_state_primaries.(PRIMARY_str{i});
    interpolators_primaries.(PRIMARY_str{i} ).spline = cell(1,6);  % For 6 dimensions (3 for position, 3 for velocity)
    for dim = 1:6
        interpolators_primaries.(PRIMARY_str{i} ).spline{dim} = spline(rtbp_times, rtbp_state_primaries.(PRIMARY_str{i} )(dim,:));
    end
end

% Combine interpolators for primaries and object into one struct
interpolators = struct();

% Add interpolators for primaries
for i = 1:length(PRIMARIES)
    PRIMARY = PRIMARIES{i};
    PRIMARY_str{i} = regexprep(PRIMARY, [{'\s+'}, {'-'}], '_');
    interpolators.(PRIMARY_str{i}) = interpolators_primaries.(PRIMARY_str{i});
end

% Add interpolators for object
for i = 1:length(object)
    BODY_str{i} = regexprep(BODY_str{i}, [{'\s+'}, {'-'}], '_');
    interpolators.(BODY_str{i}) = interpolators_object.(BODY_str{i});
end

%---------------------------------------%
cspice_kclear()
% --------------------------------------%

end


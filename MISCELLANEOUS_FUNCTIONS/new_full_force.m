function xdot = new_full_force(t_inertial,x_inertial)

global PRIMARIES
global BODIES
global FRAME
global OBSERVER
% Persistent variables
persistent GM;


    % All bodies considered for gravitational forces
    all_bodies = [PRIMARIES, BODIES];
    n_bodies = length(all_bodies);
    
    % Load required Kernels and ephemeris file once
    META = 'kernels_to_load.tm'; % Initialize required kernels
    cspice_furnsh(META); % Furnish kernels
    t_utc = cspice_et2utc(t_inertial, 'C', 0);
    cspice_kclear()
    
    % Preallocate accumulators for parallel loop
    
    a = 0;
    
    % Initialize GM if it's empty or has fewer elements than required
    if isempty(GM) || length(GM) < n_bodies
        GM = zeros(1, n_bodies);
    end

    % Parallel computation of accelerations and variational terms
    for i = 1:n_bodies
        BODY = all_bodies{i};
        
        % Check if GM is already computed for this body
        
    if GM(i) == 0
        META = 'kernels_to_load.tm'; % Initialize required kernels
        cspice_furnsh(META); % Furnish kernels
        GM(i) = cspice_bodvrd(BODY, 'GM', 1); % [km^3/s^2]
        cspice_kclear()
    end
    GM_BODY = GM(i);
    META = 'kernels_to_load.tm'; %initialize required kernels   
    cspice_furnsh(META); %furnish kernels
    [state_b, ~] = cspice_spkpos(BODY, t_inertial, 'ECLIPJ2000', 'None', OBSERVER);
    cspice_kclear()
    px = state_b(1);
    py = state_b(2);
    pz = state_b(3);

    x = x_inertial(1);
    y = x_inertial(2);
    z = x_inertial(3);
    
    r_sb = sqrt(((px - x)^2 + (py - y)^2 + (pz - z)^2));
    r_sb2 = r_sb^2;
    r_sb3 = r_sb^3; % term to speedup
    r_sb5 = r_sb^5; % term to speedup

    % Acceleration due to the bodies gravity
    a = a - GM(i) / norm( x_inertial(1:3) - state_b(1:3) )^3 * ( x_inertial(1:3) - state_b(1:3) ); % [km/s^2] 
end


% only acceleration
xdot = [x_inertial(4:6); a];






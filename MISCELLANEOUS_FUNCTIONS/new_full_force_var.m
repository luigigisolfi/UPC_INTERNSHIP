function xdot = new_full_force_var(t_inertial, x_inertial)
    % Compute the state transition matrix (STM) and the variational equations
    % for an n-body problem where the spacecraft is influenced by the gravity
    % of multiple bodies in the solar system.
    
    % Global variables
    global PRIMARIES
    global BODIES
    global OBSERVER
    
    % state_bersistent variables
    persistent GM;
    
    % All bodies considered for gravitational forces
    all_bodies = [PRIMARIES, BODIES];
    n_bodies = length(all_bodies);
    a = 0; % [km/s^2]

    
    % Load required Kernels and ephemeris file once
    META = 'kernels_to_load.tm'; % Initialize required kernels
    cspice_furnsh(META); % Furnish kernels
    %t_utc = cspice_et2utc(t_inertial, 'C', 0);
    
    % state preallocate accumulators for parallel loop
    jacobianLowerLeft = zeros(3,3);
    
    % Initialize GM if it's empty or has fewer elements than required
    if isempty(GM) || length(GM) < n_bodies
        GM = zeros(1, n_bodies);
    end
    
    % state_barallel computation of accelerations and variational terms
    for i = 1:n_bodies
        BODY = all_bodies{i};
        
        % Check if GM is already computed for this body
        if GM(i) == 0
            GM(i) = cspice_bodvrd(BODY, 'GM', 1); % [km^3/s^2]
        end

        GM_BODY = GM(i);
        three_GM = 3*GM_BODY;
        
        % Get the state of the body
        [state_b, ~] = cspice_spkpos(BODY, t_inertial, 'ECLIPJ2000', 'None', OBSERVER);
        px = state_b(1);
        py = state_b(2);
        pz = state_b(3);

        x = x_inertial(1);
        y = x_inertial(2);
        z = x_inertial(3);

        r_sb = sqrt(((px - x)^2 + (py - y)^2 + (pz - z)^2));
        r_sb3 = r_sb^3; % term to speedup
        r_sb5 = r_sb^5; % term to speedup
        
        % Acceleration due to the bodies gravity
        a = a - GM(i) / norm( x_inertial(1:3) - state_b(1:3) )^3 * ( x_inertial(1:3) - state_b(1:3) ); % [km/s^2] 
        
        % Uxx--Uzz in lower left corner
        jacobianLowerLeft(1,1) = jacobianLowerLeft(1,1) + (three_GM*(px - x)*(px - x))/(r_sb5) - GM(i)/r_sb3;
        jacobianLowerLeft(1,2) = jacobianLowerLeft(1,2) + (three_GM*(py - y)*(px - x))/(r_sb5);
        jacobianLowerLeft(1,3) = jacobianLowerLeft(1,3) + (three_GM*(pz - z)*(px - x))/(r_sb5);
        jacobianLowerLeft(2,1) = jacobianLowerLeft(1,2);
        jacobianLowerLeft(2,2) = jacobianLowerLeft(2,2) + (three_GM*(py - y)*(py - y))/(r_sb5) - GM(i)/r_sb3;
        jacobianLowerLeft(2,3) = jacobianLowerLeft(2,3) + (three_GM*(pz - z)*(py - y))/(r_sb5);
        jacobianLowerLeft(3,1) = jacobianLowerLeft(1,3);
        jacobianLowerLeft(3,2) = jacobianLowerLeft(2,3);
        jacobianLowerLeft(3,3) = jacobianLowerLeft(3,3) + (three_GM*(pz - z)*(pz - z))/(r_sb5) - GM(i)/r_sb3;
    end

    jacobian = [zeros(3),          eye(3);
                jacobianLowerLeft, zeros(3)];
    STM = jacobian * reshape(x_inertial(7:42),6,6);

    xdot = [x_inertial(4:6); a; reshape(STM,36,1)];
    cspice_kclear()
    
    
end

function xdot = new_full_force_var_vectorized(t_inertial, x_inertial)

    % Compute the state transition matrix (STM) and the variational equations
    % for an n-body problem where the spacecraft is influenced by the gravity
    % of multiple bodies in the solar system.

    % Global variables
    global PRIMARIES BODIES OBSERVER FRAME

    % Persistent variables
    persistent GM;

    % Reshape y to separate the individual trajectories

    %x_inertial = reshape(x_inertial, [], N);
    
    % All bodies considered for gravitational forces
    all_bodies = [PRIMARIES, BODIES];
    n_bodies = length(all_bodies);
    
    % clear cache and load all kernels at the first call
    persistent flagKernelLoaded;
    if isempty(flagKernelLoaded)
        cspice_kclear;
        % load all kernels from input
    % Load required Kernels and ephemeris file once
        META = './KERNELS/kernels_to_load.tm'; % Initialize required kernels
        cspice_furnsh(META); % Furnish kernels
        flagKernelLoaded = 1;
    end
    

    % Initialize GM if it's empty or has fewer elements than required
    if isempty(GM) || length(GM) < n_bodies
        GM = zeros(1, n_bodies);
    end

    % Preallocate derivative matrix
    xdot = zeros(42,1);
    % Compute accelerations and variational terms for each trajectory
        % Extract the current state for this trajectory
    
    x_current = x_inertial;

    % Initialize accumulators
    a = zeros(3, 1); % [km/s^2]
    jacobianLowerLeft = zeros(3, 3);

    x = x_current(1);
    y = x_current(2);
    z = x_current(3);
    for i = 1:n_bodies
        BODY = all_bodies{i};

        % Check if GM is already computed for this body
        if GM(i) == 0
            GM(i) = cspice_bodvrd(BODY, 'GM', 1); % [km^3/s^2]
        end

        GM_BODY = GM(i);
        three_GM = 3 * GM_BODY;

        % Get the state of the body
        [state_b, ~] = cspice_spkpos(BODY, t_inertial, FRAME, 'None', OBSERVER);
        px = state_b(1);
        py = state_b(2);
        pz = state_b(3);

        r_sb3 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(3/2); % term to speedup
        r_sb5 = ((px - x)^2 + (py - y)^2 + (pz - z)^2)^(5/2); % term to speedup

        % Acceleration due to the body's gravity
        if strcmpi(BODY,OBSERVER) 
            % center body
            a = a - GM(i) / norm(x_current(1:3) - state_b(1:3) )^3 * (x_current(1:3) - state_b(1:3) ); % [km/s^2] 
        else
            % perturbation body
            a = a - GM(i) / norm(x_current(1:3) - state_b(1:3) )^3 * (x_current(1:3) - state_b(1:3) )...
                - GM(i) / norm( state_b(1:3) )^3 * state_b(1:3); % [km/s^2]
        end
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

    jacobian = [zeros(3), eye(3);
                jacobianLowerLeft, zeros(3)];

    STM = jacobian * reshape(x_current(7:42), 6, 6);
    % Store the derivatives for this trajectory
    xdot = [x_current(4:6); a; reshape(STM, 36, 1)];
    xdot = xdot(:);
    end
    % Flatten xdot back into a column vector

    %xdot = reshape(xdot,);
end

function [mu, body_1, body_2] = get_mu(PRIMARIES)

global mu

% ----------------LOADING KERNELS --------------------------------------%
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------------------------------------------%
% Define the bodies
body_1 = PRIMARIES{1};
body_2 = PRIMARIES{end};

% Get the gravitational parameters (GM) for both bodies
GM_1 = cspice_bodvrd(body_1, 'GM', 1);
GM_2 = cspice_bodvrd(body_2, 'GM', 1);

% Determine the more massive body
if GM_1 > GM_2
    % body_1 is already the more massive body
else
    % Swap body_1 and body_2
    temp_body = body_1;
    body_1 = body_2;
    body_2 = temp_body;
    
    temp_GM = GM_1;
    GM_1 = GM_2;
    GM_2 = temp_GM;
end

% Compute the gravitational parameter mu
mu = GM_2 / (GM_1 + GM_2);

% Unload the SPICE kernel
cspice_kclear;
    
end


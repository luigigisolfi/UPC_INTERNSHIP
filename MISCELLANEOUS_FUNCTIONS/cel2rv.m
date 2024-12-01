function [r, v] = cel2rv(orbel, GM)
    %------------------------------------------------------------------
    % It computes the position r and velocity v from the classical 
    % orbital elements in the considered reference frame and units.
    % GM is the mu mass parameter [L^3/T^2] consistent with the r [L], 
    % v [L/T] output units.
    % Input:
    % orbel(1): parameter p of the orbit [L].
    % orbel(2): eccentricity
    % orbel(3): inclination (rad)
    % orbel(4): argument of the ascending node (rad)
    % orbel(5): argument of perigee (rad)
    % orbel(6): true anomaly (rad)
    %
    % Output: r, v position and velocity.
    %------------------------------------------------------------------

    % Extracting the orbital elements
    p = orbel(1);
    e = orbel(2);
    i = orbel(3);
    Omega = orbel(4);
    omega = orbel(5);
    theta = orbel(6);

    % Computing the perifocal coordinates
    costruea = cos(theta);
    sintruea = sin(theta);
    rm = p / (1 + e * costruea);
    rperif = [rm * costruea; rm * sintruea; 0];
    vperif = sqrt(GM / p) * [-sintruea; e + costruea; 0];

    % Computing the rotation matrix
    rotate = rotz(Omega * 180 / pi) * rotx(i * 180 / pi) * rotz(omega * 180 / pi);

    % Computing the position and velocity in the inertial frame
    r = rotate * rperif;
    v = rotate * vperif;
end

% Auxiliary functions for rotation matrices
function R = rotx(angle)
    % Rotation matrix around x-axis
    c = cosd(angle);
    s = sind(angle);
    R = [
        1 0 0;
        0 c -s;
        0 s c;
    ];
end

function R = rotz(angle)
    % Rotation matrix around z-axis
    c = cosd(angle);
    s = sind(angle);
    R = [
        c -s 0;
        s c 0;
        0 0 1;
    ];
end

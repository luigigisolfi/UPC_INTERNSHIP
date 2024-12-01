function orbel = rv2cel(r,v,GM)
    %------------------------------------------------------------------
    % It computes the classical orbital elements from the position r and 
    % velocity v vectors in the considered reference frame and units.
    % GM is the mu mass parameter [L^3/T^2] consistent with the r [L], 
    % v [L/T] units.
    % Input:
    % inertial state (r,v) and system's gravitational parameter (GM)
    % Output:
    % orbel(1): parameter p of the orbit [L].
    % orbel(2): eccentricity
    % orbel(3): inclination (rad)
    % orbel(4): argument of the ascending node (rad)
    % orbel(5): argument of perigee (rad)
    % orbel(6): true anomaly (rad)
    %
    % NOTES:
    % - An internal quantity (small) is considered for the case
    %   of undefined classical orbital elements. When an element 
    %   is undefined, it is set to zero, and other(s) following it are 
    %   measured according to this convention.
    % - The function stops when r and v are aligned. 
    %------------------------------------------------------------------

    small = 1.e-10;
    h = cross(r, v);
    rm = norm(r);
    vm = norm(v);
    hm = norm(h);
    
    if hm < small
        fprintf('rv2cel. rectilinear or collision trajectory. Stopping.\n');
        return;
    end
    
    n = zeros(size(h));
    n(1) = -h(2);
    n(2) = h(1);
    nm = norm(n);
    
    if nm < small
        n(1) = 1;
        n(2) = 0;
        nm = 0;
    else
        n = n / nm;
        nm = 1;
    end
    
    dotrv = dot(r, v);
    e = ((vm^2 - GM/rm) * r - dotrv * v) / GM;
    em = norm(e);
    
    if em < small
        e = n;
        em = 0;
    else
        e = e / em;
    end
    
    orbel(1) = hm^2 / GM;
    orbel(2) = em;
    aux = h(3) / hm;
    
    if abs(aux) > 1
        aux = sign(aux);
    end
    
    orbel(3) = acos(aux);
    orbel(4) = 0;
    
    if nm > small
        aux = n(1);
        if abs(aux) > 1
            aux = sign(aux);
        end
        orbel(4) = acos(aux);
    end
    
    if n(2) < 0
        orbel(4) = 2 * pi - orbel(4);
    end
    
    orbel(5) = 0.0;
    
    if em > small
        aux = dot(e, n);
        if abs(aux) > 1
            aux = sign(aux);
        end
        orbel(5) = acos(aux);
        
        if nm > small
            if e(3) < 0
                orbel(5) = 2 * pi - orbel(5);
            end
        else
            if h(3) * e(2) < 0
                orbel(5) = 2 * pi - orbel(5);
            end
        end    
    end
    
    aux = dot(e, r) / rm;
    if abs(aux) > 1
        aux = sign(aux);
    end
    
    orbel(6) = acos(aux);
    
    if em > small
        if dotrv < 0
            orbel(6) = 2 * pi - orbel(6);
        end
    else
        ha = cross(e, r);
        if dot(h, ha) < 0
            orbel(6) = 2 * pi - orbel(6);
        end
    end
end

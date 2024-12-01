function pos = get_planet_state_vector(interpolators, planet, t)

% ----------------------------------------------------------------------%
% This function gets the position of a planet at a given time by means of
% interpolation cubic splines acting on the SPICE retrieved positions. 
% ----------------------------------------------------------------------%
    
planet_str = regexprep(planet, [{'\s+'}, {'-'}], '_');
pos = zeros(1, 6);
    for dim = 1:6
        pos(dim) = ppval(interpolators.(planet_str).spline{dim}, t);
    end
end


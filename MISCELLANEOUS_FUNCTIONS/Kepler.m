function K = Kepler(E, M_e, e)
    % Kepler's equation, to be used in a Newton solver.
    K = E - e * sin(E) - M_e;
end
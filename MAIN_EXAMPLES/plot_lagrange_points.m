%----------------------------------------------------------------------------------
% This main calls the function plot_lagrange_points(mu)
% This is shown for the Earth-Moon case.
% It might be possible that Sun-Earth won't converge. In that case, 
% go to Lagrange_Points (in MISCELLANEOUS) and change the initial guesses
% for L1, L2, L3...
%----------------------------------------------------------------------------------

mu = get_mu({'MOON', 'EARTH'});
lagrange_points(mu)
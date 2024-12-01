function [xx, yy, zz] = earth_moon_sphere(varargin)
%EARTH_MOON_SPHERE Generate Earth and Moon-sized spheres.
%   [X, Y, Z] = EARTH_MOON_SPHERE(N) generates spheres representing the
%   Earth and the Moon. N determines the resolution of the spheres.
%
%   EARTH_MOON_SPHERE(N, 'unit', 'plot_option') allows specifying the unit 
%   for the axes and which object(s) to plot. Valid units are 'km', 'm', 
%   'mile', 'nm', and 'AU'. Valid plot options are 'earth', 'moon', 'both'.
%
%   Examples:
%       earth_moon_sphere
%       earth_moon_sphere(100, 'km', 'both')
%       earth_moon_sphere(50, 'mile', 'earth')
%
% Written and adapted by ChatGPT.

% Parse inputs
[cax, args, nargs] = axescheck(varargin{:});
if nargs > 3
    error('Too many input arguments. Use up to 3: resolution, units, and plot option.');
end

% Default values
n = 50; % Sphere resolution
units = 'km'; % Default unit
plot_option = 'both'; % Default plot option
for i = 1:nargs
    if isnumeric(args{i})
        n = args{i};
    elseif ischar(args{i})
        if any(strcmpi(args{i}, {'earth', 'moon', 'both'}))
            plot_option = lower(args{i});
        else
            units = args{i};
        end
    else
        error('Invalid input type. Use a number for resolution and strings for units or plot option.');
    end
end

% Scaling factors
Scale = {'moon_dist','km', 'm', 'mile', 'nm', 'au', 'ft';
         2.60e-06, 1, 1000, 0.621371192, 0.539956803, 6.68458712e-9, 3280.839895};
idx = find(strcmpi(units, Scale(1, :)), 1);
if isempty(idx)
    error('Invalid unit. Use km, m, mile, nm, au, or ft.');
end
unit_scale = Scale{2, idx};

% Earth parameters
earth_radius = 6378.1363 * unit_scale;

% Moon parameters
moon_radius = 1737.4 * unit_scale;
moon_distance = 384400 * unit_scale;

% Generate Earth sphere
[theta, phi] = meshgrid(linspace(0, 2*pi, n+1), linspace(-pi/2, pi/2, n+1));
x_earth =  0.012150582 + earth_radius * cos(phi) .* cos(theta);
y_earth = earth_radius * cos(phi) .* sin(theta);
z_earth = earth_radius * sin(phi);

% Generate Moon sphere
x_moon = moon_radius * cos(phi) .* cos(theta) - moon_distance;
y_moon = moon_radius * cos(phi) .* sin(theta);
z_moon = moon_radius * sin(phi);

% Plot Earth, Moon, or Both
if nargout == 0
    cax = newplot(cax);
    hold on;

    % Plot Earth
    if strcmp(plot_option, 'earth') || strcmp(plot_option, 'both')
        % Load Earth's topographic data if available
        if exist('topo.mat', 'file')
            load('topo.mat', 'topo');
            topo2 = [topo(:, 181:360), topo(:, 1:180)];
            surf(x_earth, y_earth, z_earth, 'FaceColor', 'texture', ...
                'CData', topo2, 'EdgeColor', 'none', 'FaceLighting', 'phong');
        else
            surf(x_earth, y_earth, z_earth, 'FaceColor', 'b', 'EdgeColor', 'none');
        end
    end

    % Plot Moon
    if strcmp(plot_option, 'moon') || strcmp(plot_option, 'both')
        surf(x_moon, y_moon, z_moon, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');
    end

    % Adjust view and labels
    axis equal;
    xlabel(['X [' units ']']);
    ylabel(['Y [' units ']']);
    zlabel(['Z [' units ']']);
    view(3);
    grid on;
    hold off;
else
    % Output coordinates
    xx = {x_earth, x_moon};
    yy = {y_earth, y_moon};
    zz = {z_earth, z_moon};
end
end

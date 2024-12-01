function [t,x] = read_orbit(filename)
% Specify the file name
% Load the data
data = load(filename);

% Assign each column to a variable
t = data(:, 1); % time
x1 = data(:, 2);
x2 = data(:, 3);
x3 = data(:, 4);
x4 = data(:, 5); %because ruilong gave me the px,py,pz momenta instead of velocities
x5 = data(:, 6);
x6 = data(:, 7);

x = [x1,x2,x3,x4,x5,x6];


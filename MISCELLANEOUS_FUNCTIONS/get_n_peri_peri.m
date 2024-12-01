function n_anomalistic = get_n_peri_peri(PRIMARIES, FRAME)

% This function computes the average value (over 100 years) of n_anomalistic for a system of
% PRIMARIES, using spice + rv2cel function

% Load MICE kernels
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels

% Define the overall time range for the search
start_time = '1950-01-01';
end_time = '2049-12-31';

% Convert to ET (Ephemeris Time)
et_start = cspice_str2et(start_time);
et_end = cspice_str2et(end_time);

% Define a coarse time step (e.g., 1 day)
time_step_coarse = 86400; % in seconds

% Generate coarse time steps
coarse_times = et_start:time_step_coarse:et_end;
num_coarse_times = length(coarse_times);

% Preallocate arrays for positions and distances
coarse_distances = zeros(1, num_coarse_times);

% Compute the Moon (or Sun's) position relative to the Earth at each coarse time step
%Computing anomalistic month duration with SPICE ephemeris

for i = 1:num_coarse_times
    [state, ~] = cspice_spkpos(PRIMARIES{1}, coarse_times(i), FRAME, 'NONE', PRIMARIES{2});
    coarse_distances(i) = norm(state);
end

% Identify local minima in the distance array (potential perigee events)
min_indices = [];
for i = 2:num_coarse_times-1
    if coarse_distances(i) < coarse_distances(i-1) && coarse_distances(i) < coarse_distances(i+1)
        min_indices(end+1) = i;
    end
end

% Preallocate for refined perigee times
refined_perigee_times = [];

% Refine each perigee event with a finer time step
for i = 1:length(min_indices)
    % Define a smaller time range around the coarse minimum
    refine_start = max(coarse_times(min_indices(i)) - time_step_coarse, et_start);
    refine_end = min(coarse_times(min_indices(i)) + time_step_coarse, et_end);
    
    % Define a finer time step (e.g., 1 hour)
    time_step_fine = 3600; % in seconds
    
    % Generate fine time steps
    fine_times = refine_start:time_step_fine:refine_end;
    num_fine_times = length(fine_times);
    
    % Preallocate arrays for fine positions and distances
    fine_distances = zeros(1, num_fine_times);
    
    % Compute the Moon (or Sun)'s position relative to the Earth at each fine time step
    for j = 1:num_fine_times
        [state, ~] = cspice_spkpos('EARTH', fine_times(j), 'J2000', 'NONE', 'MOON');
        fine_distances(j) = norm(state);
    end
    
    % Find the index of the minimum distance in the fine search
    [~, fine_min_index] = min(fine_distances);
    
    % Get the time of the minimum distance
    refined_perigee_times(end+1) = fine_times(fine_min_index); 
end

% Convert the refined perigee times back to human-readable format
n_list = zeros(length(refined_perigee_times));
refined_perigee_dates = cell(size(refined_perigee_times));
periods_list = diff(refined_perigee_times/86400);
n_list = (2*pi./mean(periods_list));
n_anomalistic = n_list/86400;
% Unload SPICE kernels

cspice_kclear()

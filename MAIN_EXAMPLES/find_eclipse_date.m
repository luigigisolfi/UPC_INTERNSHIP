% Initialize required kernels
META = './KERNELS/kernels_to_load.tm';
cspice_furnsh(META);  % Load necessary SPICE kernels

% Define the time window
start_time = '2000 Jan 01 11:58:56';
end_time = '2000 Jan 30 00:00:00';

% Convert time to ephemeris time (ET)
et_start = cspice_str2et(start_time);
et_end = cspice_str2et(end_time);
dt = 1000;  % Time step in seconds

% Initialize an empty array to store the times of detected eclipses
result = [];

% Search for occultation within the specified time window
for et = et_start : dt : et_end
    % Check for occultation of the Moon by Earth relative to the Sun
    ocltid = cspice_occult('MOON', 'ELLIPSOID', 'IAU_MOON', ...
                           'EARTH', 'ELLIPSOID', 'IAU_EARTH', ...
                           'NONE', 'SUN', et);
                       
    % We’re looking for a "total" or "annular" occultation
    if ocltid == -3
        % Store the time of eclipse in the result array
        result = [result; et];
    end
end

% Initialize an empty cell array to store unique eclipse dates
eclipse_dates = {}; 

% Loop through the results to convert ET to UTC
for i = 1:length(result)
    % Convert the ephemeris time to UTC format for each event
    utc_time = cspice_et2utc(result(i), 'C', 3); % Convert to UTC
    
    % Create a datetime object from the UTC time string
    dt = datetime(utc_time, 'InputFormat', 'yyyy MMM dd HH:mm:ss.SSS');
    
    % Extract the date part only
    date_only = datestr(dt, 'yyyy-mm-dd'); % Format the date as needed
    
    % Append the date to the array (if not already present)s
    if ~ismember(date_only, eclipse_dates)
        eclipse_dates{end + 1} = date_only; % Add the new unique date
    end
end

% Check if we found any unique eclipse dates and display them
if ~isempty(eclipse_dates)
    fprintf('Lunar eclipse occurred on the following date(s):\n');
    for i = 1:length(eclipse_dates)
        fprintf('%s\n', eclipse_dates{i}); % Display each unique date
    end
else
    fprintf('No lunar eclipse found in the specified time range.\n');
end

% Unload kernels
cspice_kclear;

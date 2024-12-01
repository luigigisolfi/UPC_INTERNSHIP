function L = get_L(PRIMARIES, FRAME)

% -----------------------------------------------------------------%
% This function computes the mean distance of two given primaries over a
% timespan of 100 years 
% -----------------------------------------------------------------%

ti = '1950 JAN 2 TDB';
tf = '2049 DEC 31 TDB';

et_start = cspice_str2et(ti);
et_end = cspice_str2et(tf);

time_step_coarse = 86400; % in seconds
coarse_times = et_start:time_step_coarse:et_end;
num_coarse_times = length(coarse_times);
for i = 1:num_coarse_times
    [state, ~] = cspice_spkpos(PRIMARIES{2}, coarse_times(i), FRAME, 'NONE', PRIMARIES{1});
    coarse_distances(i) = norm(state);
end

L = mean(coarse_distances);

end
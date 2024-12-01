function avg_e = get_avg_e(PRIMARIES, FRAME)

% This function computes the average value (over 100 years) of eccentricity for a system of
% PRIMARIES, using spice + rv2cel function

ti = '1950 JAN 2 TDB';
tf = '2049 DEC 31 TDB';

et_start = cspice_str2et(ti);
et_end = cspice_str2et(tf);

GM = get_GM_body(PRIMARIES{1});
time_step_coarse = 86400; % in seconds
coarse_times = et_start:time_step_coarse:et_end;
num_coarse_times = length(coarse_times);
e_list = zeros(num_coarse_times,1);

META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels

META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
for i = 1:num_coarse_times
    [state, ~] = cspice_spkezr(PRIMARIES{2}, coarse_times(i), FRAME, 'NONE', PRIMARIES{1});
    r = state(1:3);
    v = state(4:6);
    orbel = rv2cel(r, v, GM);
    e_list(i) = orbel(2);
end

avg_e = median(e_list);
end

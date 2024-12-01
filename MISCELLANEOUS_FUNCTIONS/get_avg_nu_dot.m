function avg_nu_dot = get_avg_nu_dot(avg_e, inertial_state_primaries, PRIMARIES, L)
global n

    GM_1 = get_GM_body(PRIMARIES{1}); %of the central body (in the case SUN-EM Barycenter, it is the SUN)
    GM_2 = get_GM_body(PRIMARIES{2});
    PRIMARY_str_2 = regexprep(PRIMARIES{2}, [{'\s+'}, {'-'}], '_');
    r_2_list = inertial_state_primaries.(PRIMARY_str_2).position; %position of PRIMARY{2}
    v_2_list = inertial_state_primaries.(PRIMARY_str_2).velocity;%velocity of PRIMARY {2}

nu_list = zeros(length(r_2_list),1);
for i = 1:length(r_2_list)
    orbel = rv2cel(r_2_list(1:3,i), v_2_list(1:3,i), GM_1);
    nu_list(i) = orbel(6);
end

nu_list(1);
nu_dot_list =sqrt(1 + avg_e*cos(nu_list))*n;
avg_nu_dot = mean(nu_dot_list);
end

function nu_eclipse = get_nu_eclipse(PRIMARIES, inertial_t0)

% Computes the true anomaly of the secondary wrt primary at a time t0
% (useful, for instance, to know true anomaly at eclipse time)

META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
inertial_state_body = cspice_spkezr(PRIMARIES{2}, inertial_t0, 'J2000', 'None', PRIMARIES{1});
r = inertial_state_body(1:3);
v = inertial_state_body(4:6);
GM = get_GM_body(PRIMARIES{1});
orbel = rv2cel(r,v,GM);
nu_eclipse = orbel(6);
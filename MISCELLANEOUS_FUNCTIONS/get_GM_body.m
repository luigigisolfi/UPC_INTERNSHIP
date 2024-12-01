function GM_body = get_GM_body(body)

META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
[val] = cspice_bodvrd(body, 'GM', 1);
%fprintf('GM of %s: %.5e (in km^3/s^2)\n', body, val);
GM_body = val;

cspice_kclear()


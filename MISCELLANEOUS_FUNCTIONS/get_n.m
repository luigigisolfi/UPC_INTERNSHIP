function n = get_n(PRIMARIES, FRAME)

% ----------------------------------------------------------------------%
% This function uses the kernel Gravity.tpc and computes the mean motion
% of a given system made of two PRIMARIES
% ----------------------------------------------------------------------%

% ----------------LOADING KERNELS --------------------------------------%
META = './KERNELS/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------------------------------------------%

% -----------------------------------------------%
% The user is allowed to input the PRIMARIES in every order 
% This is taken care of in the following, computing the max_GM
GM_1 = cspice_bodvrd(PRIMARIES{1}, 'GM', 1);
GM_2 = cspice_bodvrd(PRIMARIES{2}, 'GM', 1);
L = get_L(PRIMARIES, FRAME);
n = sqrt((GM_1 + GM_2)/L^3);
% ----------------------------------------------------------------------%
cspice_kclear()
% ----------------------------------------------------------------------%

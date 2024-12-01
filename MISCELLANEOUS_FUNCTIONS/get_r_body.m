function x_b = get_r_body(et, body, FRAME, OBSERVER)

% t is in the format 'YYYY MM DD'%
% all inputs are strings %
%body = string(body)
fprintf('EHYYYYY %s', body)
TDBFMT = 'YYYY MON DD HR:MN:SC.### (TDB) ::TDB';
%et = cspice_str2et(t);
timstr = cspice_timout(et, TDBFMT );
body;
FRAME;
OBSERVER;
%fprintf( '   Computing position of %body   at time %et               = %s\n', body, et, timstr );
[state_body, ltime_body] = cspice_spkezr(body, et, FRAME, 'NONE', OBSERVER);
x_b = state_body(1:6); %in Km;






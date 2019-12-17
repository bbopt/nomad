function [Zh,Sh,Zv,Sv] = sgtelib_server_cv

% sgtelib_server_ping;

% Remove flags
system('rm -f flag_cv_* 2>/dev/null');

% Create flag
system('touch flag_cv_transmit');

% Wait for reception flag
sgtelib_server_wait_file('flag_cv_finished');

% Read Output file
[Zh,Sh,Zv,Sv] = sgtelib_server_read_matrix('flag_cv_finished');

% Remove all flags
system('rm -f flag_cv_* 2>/dev/null');


function M = sgtelib_server_info

% sgtelib_server_ping;

% Remove flags
system('rm -f flag_info_* 2>/dev/null');

% Write infoion point
system('touch flag_info_transmit');

% Wait for reception flag
sgtelib_server_wait_file('flag_info_finished');

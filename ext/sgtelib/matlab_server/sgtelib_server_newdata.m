function sgtelib_server_newdata(X,Z)

% sgtelib_server_ping;

% Remove all flags
system('rm -f flag_new_data_* 2>/dev/null');

% Write matrices
sgtelib_server_write_matrix(X,'X','new_data_x.txt');

sgtelib_server_write_matrix(Z,'Z','new_data_z.txt');

% Create flag
system('touch flag_new_data_transmit');

% Wait for reception flag
sgtelib_server_wait_file('flag_new_data_received');

% Remove all flags
system('rm -f flag_new_data_* 2>/dev/null');

function [Z,std,ei,cdf] = sgtelib_server_predict(X)

% sgtelib_server_ping;

% Remove flags
system('rm -f flag_predict_* 2>/dev/null');

% Write prediction point
sgtelib_server_write_matrix(X,'X','flag_predict_create');

% Create flag
system('mv flag_predict_create flag_predict_transmit');

% Wait for reception flag
sgtelib_server_wait_file('flag_predict_finished');

% Read Output file
[Z,std,ei,cdf] = sgtelib_server_read_matrix('flag_predict_finished');

% Remove all flags
system('rm -f flag_predict_* 2>/dev/null');


function M = sgtelib_server_metric(metric_str)

% sgtelib_server_ping;

% Remove flags
system('rm -f flag_metric_* 2>/dev/null');

% Write metricion point
system(['echo ' metric_str ' >> flag_metric_create']);

% Create flag
system('mv flag_metric_create flag_metric_transmit');

% Wait for reception flag
sgtelib_server_wait_file('flag_metric_finished');

% Read Output file
M = importdata('flag_metric_finished');

% Remove all flags
system('rm -f flag_metric_* 2>/dev/null');

function ready = sgtelib_server_ping

system('touch flag_ping');

% Wait for reception flag

i = sgtelib_server_wait_file('flag_pong',0.1);
if i==0
    disp('Retry ping...');
    system('touch flag_ping');
    i = sgtelib_server_wait_file('flag_pong',0.5);
end

if i
    ready = importdata('flag_pong');
    system('rm -f flag_pong');
    %disp('ping ok!');
else
    disp('=====================SGTELIB_SERVER ERROR==========================');
    disp('sgtelib_server not responding');
    error('We tried to "ping" sgtelib_server, but it is off or not responding');
end



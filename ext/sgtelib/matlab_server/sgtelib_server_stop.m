function sgtelib_server_stop
disp('Kill sgtelib_server.exe');
!touch flag_quit
%!killName sgtelib_server.exe

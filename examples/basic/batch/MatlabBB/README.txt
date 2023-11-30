% WARNING: The command to run matlab in batch mode may vary with Matlab version and OS.
% The given bat (bb.bat) command has been tested with Matlab 2021b on Windows. 
%     - To make sure the command works properly with your Matlab version, start a Windows shell command terminal.
%     - Go into the directory of this example 
%     - Test the following command: bb.bat X0.txt
%     - Matlab should start in background and run the 'fun.m' function.
%     - The command should only display '4 -10 -30', that is the objective and constraints values for X0.
%     - When the blackbox Matlab command is working you can run Nomad on it.
% 
% The given shell script command (bb.sh) has been tested with Matlab 2022a on OSX and Linux
%     - To test the command, go into the directory of this example 
%     - Run the command: ./bb.sh X0.txt
%     - Matlab should start in background and run the 'fun.m' function.
%     - The command should only display '4 -10 -30', that is the objective and constraints values for X0.

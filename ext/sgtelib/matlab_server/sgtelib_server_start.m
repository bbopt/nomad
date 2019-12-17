function sgtelib_server_start(model,keepopen)

sgtelib_server_stop;
if nargin==1
    keepopen=false;
end

disp('Start sgtelib.exe in server mode.');


% Selection of the terminal software
TERMINAL_SOFTWARE_LIST = { 'gnome-terminal' 'lxterm' 'uxterm' 'xterm' 'konsole'};
for i=1:length(TERMINAL_SOFTWARE_LIST)
    termprog = TERMINAL_SOFTWARE_LIST{i};
    if ~system(['which ' termprog]);
        break;
    else
        disp(['Could not find terminal software: ' termprog]);
        termprog = 'bg';
    end
end

termprog = 'bg'

disp(['Selected terminal software: ' termprog]);


% Option to start sgtelib in gdb
gdboption = ' ';
% Verbose option of sgtelib.
verboseoption = ' ';
% Option of keep open
if keepopen
    if ~system('which gdb')
        gdboption = ' gdb -q -ex run --args ';
    end
    verboseoption = ' -verbose ';
end
% command to start sgtelib.
sgtelibcmd = [' sgtelib.exe -server -model ' model verboseoption];


% Reset ld_library_path
old_ld_library_path = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH','.');

switch termprog
    case 'gnome-terminal'
        termoption = ' -t sgtelib_server --hide-menubar -e ';
        % Command to run after the end of sgtelib
        if keepopen
            endoption = ' ; exec /bin/bash -i ';
        else
            endoption = ' ';
        end
        command = [termprog termoption '"/bin/bash -c '' ' gdboption sgtelibcmd endoption '''" &'];
    case {'xterm','lxterm','uxterm'}
        if keepopen
            termoption = ' -hold -e ';
        else
            termoption = ' -e ';
        end
        command = [termprog termoption gdboption sgtelibcmd ' &'];
    case 'konsole'
        if keepopen
            termoption = ' --hold -e ';
        else
            termoption = ' -e ';
        end
        command = [termprog termoption gdboption sgtelibcmd ' &'];
    case 'bg'
        command = [sgtelibcmd ' &'];
end

disp(command)
[status,response] = system(command);
if keepopen
    pause(1);
end
if status || ~isempty(response)
    disp(command);
    disp(status)
    disp(response)
end
pause(1);

% Old LD_LIBRARY_PATH
setenv('LD_LIBRARY_PATH',old_ld_library_path);


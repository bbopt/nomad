close all
clear all
clc

% Read the data in help_data.tex and convert it to build the file sgtelib_help.cpp


% Commands to look for
COMMANDS = {'\itemname','\itemkw','\iteminfo','\section','\subsection','\helpdivider'};
% Total number of commands
NC = length(COMMANDS);
% Number of commands that shoult not be ignored (the 3 first commands)
NC_OK = 3;

% Cell containint all the help data
DATA = cell(0,NC_OK);

findcmd = false(NC,1);
str = '';
nc = 1;
cj = 0;

% Open the tex file
file = 'help_data.tex';
fid = fopen(file);
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end

    disp(line);
    % Clean the line
    line = cleanSpaces(line);

    if isempty(line) || line(1)=='%'
        continue;
    end

    % Find which command(s) are in the line
    for i=1:NC
        findcmd(i) = ~isempty(strfind(line,COMMANDS{i}));
    end
    % Error if more than one command.
    if sum(findcmd)>1
        disp(line);
        error('Several commands in this line!');
    end
    % If there is exactly one command
    if sum(findcmd)==1
        % Update the index of the command.
        disp('===================');
        nc = find(findcmd)
        if nc==1
            cj = cj+1;
        end
        cj
        % If the command is acceptable
        if nc<=NC_OK
            % Open new string in DATA
            DATA{cj,nc} = '';
        end
        
        % Remove the command from the line
        line = strrep(line,COMMANDS{nc},'');
    end

    if nc<=NC_OK
        DATA{cj,nc} = [DATA{cj,nc} ' ' line];
    end

end
fclose(fid);

% List of strings to replace in the data
REPLACES = { 
    '\ref{',''
    '\spe{',''
    '\newline','\n'
    '\linebreak','\n'
    '\\','\n'
    '\;',' '
    '\,',' '
    '\example{','     '
    '$<$','<'
    '$>$','>'
    '$\ge$','>='
    '$\le$','<='
    '\item','\n *'
    '\begin{itemize}',''
    '\end{itemize}','\n'
    '\quote{{\tt ',''
    '{\tt ',''
    '}',''
    '"','\"'
    '\myUnderline{',''
    '\smalltitle{','\n'
    '\textbf{',''
    '\textbf{',''
    '\_','_'
    '{',''
    '}',''
    '$',''
    '\n','SPECIALKWNEWLINE'
    '\"','SPECIALKWQUOTE'
    '\',''
    'SPECIALKWNEWLINE','\n'
    'SPECIALKWQUOTE','\"'
    };


% Clean the data
for cj=1:size(DATA,1)
    for nc=1:NC_OK
        % Clean string
        str = cleanSpaces(DATA{cj,nc});
        % Remove brackets.
        if str(1) == '{'
            str(1) = [];
        else
            disp(str);
            error('Does not start with curbe bracket');
        end
        if str(end) == '}'
            str(end) = [];
        else
            disp(str);
            error('Does not end with curbe bracket');
        end

        % Clean string again
        str = cleanSpaces(str);
        % Other cleaning
        for k=1:size(REPLACES,1)
            str = strrep(str,REPLACES{k,1},REPLACES{k,2});
        end
        while length(str)>=2 && strcmp(str(1:2),'\n')
            str(1:2) = [];
        end
        while length(str)>=2 && strcmp(str(end-1:end),'\n')
            str(end-1:end) = [];
        end
        str = strrep(str,'\n',['\n"' char(10) '"'])
        % Other modifications
        DATA{cj,nc} = str;
        if nc==1
            DATA{cj,nc} = upper(DATA{cj,nc});
        end
    end
end

NL = size(DATA,1);
for i=1:NL
    for j=1:NC_OK
        disp([num2str(i) ',' num2str(j) '==' DATA{i,j} '==']);
    end
end


% Write cpp file.
file = '../src/sgtelib_help.cpp';

copyfile('help_header.cpp',file);

fid = fopen(file,'a');
newline = char(10);
mywrite = @(line) fwrite(fid,[line newline]);

mywrite('//================================');
mywrite('//  Get dimension of HELP_DATA');
mywrite('//================================');
mywrite('int SGTELIB::dim_help_data (void){');
mywrite(['  return ' num2str(NL) ';']);
mywrite('}//');
mywrite('');
mywrite('//================================');
mywrite('//  Construct the help data');
mywrite('//================================');
mywrite('std::string ** SGTELIB::get_help_data (void){');
mywrite = @(line) fwrite(fid,['  ' line newline]);
mywrite('int i;');
mywrite(['const int NL = ' num2str(NL) ';']);
mywrite(['const int NC = ' num2str(NC_OK) ';']);
mywrite('std::string ** HELP_DATA = new std::string * [NL];');
mywrite('for (i = 0 ; i<NL ; i++) HELP_DATA[i] = new std::string [NC];');
mywrite('i = 0;');

for i=1:NL
    mywrite(['//================================']);
    mywrite(['//      ' DATA{i,1}]);
    mywrite(['//================================']);   
    for j=1:NC_OK
        mywrite(['HELP_DATA[i][' num2str(j-1) '] = "' DATA{i,j} '";']);

    end
    mywrite(['i++;']);
end

mywrite(['//================================']);
mywrite('return HELP_DATA;');

mywrite = @(line) fwrite(fid,[line newline]);
mywrite('}//');

fclose(fid);
edit(file);



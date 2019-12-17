function [varargout] = sgtelib_server_read_matrix(file)


% Get matrices names from the file
NAMES = cell(0,0);
fid = fopen(file);
line = fgetl(fid);
while ischar(line)
    if ~isempty(line)
        if any(line=='=')
            i = find(line=='=');
            line = line(1:i-1);
            line=strrep(line,' ','');
            NAMES{end+1} = line;
        end
    end
    line = fgetl(fid);
end

N = length(NAMES);

% Check that the number of output is smaller than the number
% of matrices in the file
if N<nargout
    disp(['File name: ' file]);
    disp(['Nb of matrices in the file: ' num2str(N)]);
    disp(['Nb of matrices required in output: ' num2str(nargout)]);
    error('The file does not contain enough matrices');
end

% Convert the file in .m file and evaluate it.
mfile = ['matlab_' file num2str(rand,12)];
mfile = strrep(mfile,'.','');
mfile = strrep(mfile,'-','');
copyfile(file,[mfile '.m']);
eval(mfile);
delete([mfile '.m']);

% Affect the matrices in varargout
for i=1:N
    varargout{i} = eval(NAMES{i});
end


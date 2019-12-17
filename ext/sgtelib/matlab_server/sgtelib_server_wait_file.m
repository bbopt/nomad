function i = sgtelib_server_wait_file(name,wait_tmax)

if ~iscell(name)
    name = {name};
end

if nargin==1
    wait_tmax = 1000; 
end 
wait_dt = 0.001;
wait_t = 0;

while wait_t < wait_tmax
  pause(wait_dt);
  wait_t = wait_t+wait_dt;

  % Check if there is a new input file
  for i=1:length(name)
      if exist(name{i},'file')
          return;
      end
  end
end

i = 0;
if length(name)==1
    disp(['sgtelib_server_wait_file: file "' name{1} '" not found within time limit']);
else
    disp(['sgtelib_server_wait_file: file not found within time limit']);
end

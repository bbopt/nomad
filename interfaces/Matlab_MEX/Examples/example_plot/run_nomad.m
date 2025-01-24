% CRESCENT 5

clear fun; % Clear persistent index

x0 = [0 0 0 0 0 ]';
lb = [-6 -6 -6 -6 -6 ]';
ub = [6 6 6 6 6]';
params = struct('display_degree','2','max_bb_eval','100','direction_type','ortho 2n','bb_output_type','OBJ PB PB'); 

% Create the plot
figure;
hold on;

% Start optimization
[x,fval] = nomadOpt(@fun,x0,lb,ub,params);


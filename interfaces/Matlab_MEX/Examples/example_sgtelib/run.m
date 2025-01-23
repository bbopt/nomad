% Rosenbrock
%clc

x0 = [0;0];
params = struct('display_degree','2','max_eval','100','bb_output_type','OBJ','param_file','param.txt');

[xr,fval,ef,iter] = nomadOpt(@bb,x0,[-10;-10],[10;10],params);



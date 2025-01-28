%% MINLP 1 [fval = -5]
%clc
fun = @(x) [ (x(1) - 5)^2 + x(2)^2 - 25; 
              x(1)^2 - x(2) + 0.5 ]';
x0 = [10; 10];
params = struct('display_degree','2','bb_input_type','*I','max_eval','100','direction_type','ortho 2n','bb_output_type','OBJ PB');

[xr,fval,ef,iter] = nomadOpt(fun,x0,[],[],params);



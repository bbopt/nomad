% Example for user callback to stop optimization 

fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [-2 1]';
lb = [-Inf;-1.5];
ub = [100;100];


params = struct('initial_mesh_size','* 10','MAX_BB_EVAL','100','DISPLAY_ALL_EVAL','yes','DISPLAY_STATS','BBE ( SOL ) BBO'); 

% Start optimization
[x,fval,hinf,exit_status,nfeval] = nomadOpt(fun,x0,lb,ub,params,@cbFun);


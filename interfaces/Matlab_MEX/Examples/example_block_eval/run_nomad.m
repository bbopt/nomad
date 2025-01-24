% CRESCENT 5

% Parallel pool must be started separately (for multiple runs) or in the
% script (single run)
%parpool(2);

x0 = [0 0 0 0 0 ];
lb = [-6 -6 -6 -6 -6 ];
ub = [6 6 6 6 6];

params = struct('display_degree','2','max_bb_eval','100','direction_type','ortho 2n','bb_output_type','OBJ PB PB','bb_max_block_size','4'); 

% Start optimization

% Without callback
% [x,fval,hinf,exit_status,nfeval] = nomadOpt(@fun,x0,lb,ub,params);

% With callback
[x,fval,hinf,exit_status,nfeval] = nomadOpt(@fun,x0,lb,ub,params,@cbFun);

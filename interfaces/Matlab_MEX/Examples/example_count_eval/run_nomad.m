% This example illustrates the case where blackbox evaluations are not 
% counted when a hidden constraint is not satisfied. See bb.m.

clc
x0 = [1;1];

params = struct('display_all_eval','yes','display_degree','2','bb_output_type','CNT_EVAL OBJ', 'MAX_BB_EVAL','30','display_stats','bbe bbo');
[x,fval] = nomadOpt(@bb,x0,[-10;-10],[10;10], params);

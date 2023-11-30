%  F. Kursawe, “A variant of evolution strategies for vector optimization,” in PPSN I, Vol 496 Lect Notes in Comput Sc. Springer-Verlag, 1991, pp. 193–197.

x0 = [0 0 0 ];
lb = [-5 -5 -5 ];
ub = [5 5 5 ];

params = struct('display_degree','2','max_eval','1000','direction_type','ortho 2n','bb_output_type','OBJ OBJ','DMULTIMADS_OPTIMIZATION','true','SOLUTION_FILE','sol.txt');

[x,fval,hinf,exit_status,nfeval] = nomadOpt(@bb,x0,lb,ub,params);



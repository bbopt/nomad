x0 = [0 100 1 100 ]';
lb = [-100 0.0 0 0.0  ]';
ub = [100 10000.0 100 10000]';
opts = nomadset('display_degree',2,'display_all_eval',1,'history_file','history.txt','max_bb_eval',500,'bb_output_type','EB EB OBJ','f_target',0,'bb_input_type','[C R C R]','neighbors_mat',@neighbors); 

% Start optimization
[x,fval] = nomad(@fun,x0,lb,ub,opts);

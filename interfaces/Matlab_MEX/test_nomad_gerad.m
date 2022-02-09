% GERAD VERSION TESTING


%% PROBLEM 1
fun = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
x0 = [-2 1]';
lb = [-Inf;-1.5];
ub = [100;100];

params = struct('initial_mesh_size','* 10','MAX_BB_EVAL','100'); 

% Start optimization
[x,fval,hinf,exit_status,nfeval] = nomadOpt(fun,x0,lb,ub,params);


%Uncomment the following problems for further testing


% %% PROBLEM 2
% % Blackbox Function
% %clc
% bb = @(x) [29.4*x(1) + 18*x(2);
%            -(x(1) - 0.2458*x(1)^2/x(2)) + 6];
% % Bounds      
% lb = [0;1e-5];
% ub = [115.8;30];
% % Starting Guess
% x0 = [0;1e-5];
% % Options
% params = struct('display_degree','2','bb_output_type','OBJ EB','max_bb_eval','50');
% 
% [x,fval,hinf,exit_status,nfeval] =  nomadOpt(bb,x0,lb,ub,params);
% 
% 
% 
% 
% %% PROBLEM 3
% %clc
% % Blackbox Function
% bb = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% % Starting Guess
% x0 = [0 0]';
% params = struct('display_degree','2', 'direction_type','ortho n+1 neg', 'max_bb_eval','200');
% % Solve
% [x,fval,hinf,exit_status,nfeval] = nomadOpt(bb,x0,[],[],params);
% 
% PROBLEM 4 [fval = -2.5]
% %clc
% fun = @(x)   [-x(1) - x(2) - x(3);
%               (x(2) - 1./2.)*(x(2) - 1./2.) + (x(3) - 1./2.)*(x(3) - 1./2.) - 1/4;
%                 x(1) - x(2);
%                 x(1) + x(3) + x(4) - 2];          
% ub = [1;10;10;5];
% lb = [0;0;0;0];
% x0 = [0;0;0;0];
% 
% params = struct('display_degree','2', 'bb_output_type','OBJ PB PB PB','bb_input_type','(B R R I)', 'max_bb_eval','50', 'eval_queue_sort','dir_last_success');
% 
% [x,fval,hinf,exit_status,nfeval] = nomadOpt(fun,x0,lb,ub,params);
% 
% 
% %% Rosenbrock [x = 1,1, fval = 0]
% % Blackbox Function
% %clc
% bb = @(x) (1-x(1))^2 + 100 *(x(2)-x(1)^2)^2;
% % Starting Guess
% x0 = [0.5 0]';
% % Solve
% [x,fval,ef,iter] = nomadOpt(bb,x0)
% 
%% St_e01 [x = 6,0.6667, fval = -6.6667]
% % Blackbox Function
% %clc
% bb = @(x) [-x(1) - x(2);
%            x(1)*x(2) - 4];
% % Bounds
% lb = [0;0];
% ub = [6;4];
% % Starting Guess
% x0 = [1,1]';
% % Options
% params = struct('display_degree','2','bb_output_type','OBJ PB');
% % Solve
% [x,fval,ef,iter] = nomadOpt(bb,x0,lb,ub,params)
% 
% %% St_e08 [x = 0.1294, 0.4830, fval = 0.7418]
% %clc
% % Blackbox Function
% bb = @(x) [2*x(1)+x(2);
%            -16*x(1)*x(2) + 1; 
%            (-4*x(1)^2) - 4*x(2)^2 + 1];
% % Bounds       
% lb = [0;0];
% ub = [1;1];
% % Starting Guess
% x0 = [0;0];
% % Options
% params = struct('display_degree','2','bb_output_type','OBJ PB PB');
% 
% [x,fval,ef,iter] = nomadOpt(bb,x0,lb,ub,params);
% 
%% Wolfram problem (multiple global minima fval = 0)
% %clc
% % Fitting Function
% obj = @(x) [x(1) - sin(2*x(1) + 3*x(2)) - cos(3*x(1) - 5*x(2));
%           x(2) - sin(x(1) - 2*x(2)) + cos(x(1) + 3*x(2))];      
% % Blackbox Function      
% bb = @(x) norm(obj(x));
% 
% % Bounds
% lb = [-4;-4];
% ub = [4;4];
% % Starting Guess
% x0 = [0;-3];
% % Options
% params = struct('display_degree','2','vns_mads_search','true');
% 
% [x,fval,ef,iter] = nomadOpt(bb,x0,lb,ub,params);
% 
% %% MINLP 1 [fval = -5]
% clc
% fun = @(x) [ (x(1) - 5)^2 + x(2)^2 - 25;
%               x(1)^2 - x(2) + 0.5 ];
% x0 = [0;0];
% params = struct('display_degree','2','bb_input_type','(I I)','bb_output_type','OBJ PB','max_bb_eval','20');
% [xr,fval,ef,iter] = nomadOpt(fun,x0,[],[],params)
% % 

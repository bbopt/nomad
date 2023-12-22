%% TEST SURROGATE
clc
x0 = [-10;-10];

% Let's disable quad model search to show the impact of trial point sorting more clearly

% Default is to use quad model for sorting trial points before true objective evaluations.
% params = struct('display_degree','2','bb_output_type','OBJ','MAX_BB_EVAL','1000','QUAD_MODEL_SEARCH','false');

% Use "surrogate" for sorting trial points before true objective evaluations
params = struct('display_degree','2','bb_output_type','OBJ','MAX_BB_EVAL','1000','QUAD_MODEL_SEARCH','false','EVAL_QUEUE_SORT','SURROGATE');
[x,fval] = nomadOpt(@bbsur,x0,[-10;-10],[10;10],params);

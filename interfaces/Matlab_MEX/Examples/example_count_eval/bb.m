function [eval ] = bb(x)

% IMPORTANT
% This example illustrates the case where blackbox evaluations are not 
% counted when a hidden constraint is not satisfied. 

% Let's say this is a hidden constraint (i.e. not in BB_OUTPUT_TYPE)
g = 1-x(1);

if ( g > 0 )
    eval = [0 , Inf] ; % The evaluation is not counted
else
    eval =[1 , 10*(x(2)-x(1)^2)];
end


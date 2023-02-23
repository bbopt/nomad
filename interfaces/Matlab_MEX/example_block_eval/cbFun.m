function stop = cbFun(count,out,x)

nb_points=size(x,1); % Number of evaluations
n=size(x,2); % Dimension of each point

nbOutValue = size(out,1); % Number of evaluations
n_out=size(out,2); % Number of constraints + objective

if ( n ~= 5)
    error('Black box function must be provided with points of dimension 5 (=number of columns)');
end

% Take the best of all objective function values (BB_OUTPUT_TYPE defines
% the type and order of output value)
[best_obj, index]=min(out(:,1));
% Feasibility of the best objective function value
c1_forBest = out(index,2);
c2_forBest = out(index,3);

stop = false;
if (count > 10 && best_obj < -2.0 && c1_forBest <=0 && c2_forBest <= 0)
    stop = true;
end
return

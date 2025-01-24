function [ fout ] = fun( X )

% Share index between different fun call
persistent index
if isempty(index)
    index = 0;
end

n=size(X,2);

if ( n ~= 5)
    error('Black box function must be provided with points of dimension 5');
end

% CRESCENT 5 
%%%%%%%%%%%%%%

% objective
obj=X(5);

% constraints
c1=(sum(X(:)-1).^2)-25;
c2=25-sum((X(:)+1).^2);

% Outputs
fout=[obj , c1 , c2];



% update the index 
index = index +1;

% Plot only feasible point
if ( c1 < 0 && c2 < 0)
 plot (index,obj,'bo');
 drawnow;
end
function [ fout ] = fun( X )

% CRESCENT 5 --- Block evals version

nb_points=size(X,1);
n=size(X,2);

if ( n ~= 5)
    error('Black box function must be provided with points of dimension 5 (=number of columns)');
end

% objective
obj=X(:,5);
% constraints
for m=1:nb_points
    c1(m)=sum((X(m,:)-1).^2)-25;
    c2(m)=25-sum((X(m,:)+1).^2);
end

% Unlike in the standard nomad version, a point cannot be
% rejected in the matlab version.
% If some points cannot be evaluated, objective and constraints must be set
% to Inf, otherwise the nomad evaluator for block of evaluations cannot
% continue. 
fout=[obj , c1' , c2'];


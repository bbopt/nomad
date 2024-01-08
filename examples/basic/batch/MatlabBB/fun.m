
% CRESCENT 5

X = load('inputMatlab.txt');

n=size(X,2);

if ( n ~= 5)
    error('Black box function must be provided with points of dimension 5 (=number of columns)');
end

% objective
obj=X(5);
% constraints
c1=sum((X(:)-1).^2)-25;
c2=25-sum((X(:)+1).^2);

disp([obj c1 c2]);


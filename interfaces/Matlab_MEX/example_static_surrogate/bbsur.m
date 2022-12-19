function eval = bbsur(x,sur)

if (nargin==1)
    sur=false;
end

% Almost same function (just some scaling) for surrogate and truth (see the impact on performance with or without using surrogate)
if ( sur == true )
    eval = [10*(x(2)-x(1)^2)^2+0.1*(1-x(1))^2];
else
    eval = [100*(x(2)-x(1)^2)^2+(1-x(1))^2];
end

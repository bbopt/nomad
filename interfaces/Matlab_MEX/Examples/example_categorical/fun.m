function [ out ] = fun( X )

t0=X(1);
v0=X(2);
t1=X(3);
v1=X(4);

v=zeros(3,1);
vmin = 10000;

if ( t0 < 0 || t1 < 0 )
    out(1)=Inf;
    out(2)=Inf;
    out(3)=Inf;
    return;
end

v(t0+1) = v0;
v(t1+1) = v1;

if ( v0 < vmin )
    vmin = v0;
end
if ( v1 < vmin )
    vmin = v1;
end

% constraints (budget and each asset is considered with at least 1$):
vt = v(1) + v(2) + v(3); 
h  = vt - 10000;
g1 = h;
g2 = 1-vmin;

if ( h <= 0 && vmin >= 1 )
%% compute the risk and revenue:
    vt2  = vt*vt;
    rev  = v(1) * 0.0891 + v(2) * 0.2137 + v(3) * 0.2346;
    risk = 0.01 * (v(1)/vt)*(v(1)/vt) + ...
           0.05 * (v(2)/vt)*(v(2)/vt) + ...
           0.09 * (v(3)/vt)*(v(3)/vt) + ...
           0.02 * (v(1)*v(2)/vt2)     + ...
           0.02 * (v(1)*v(3)/vt2)     + ...
           0.10 * (v(2)*v(3)/vt2);

    % the objective is taken as a scaled distance
    % between (risk,revenue) and (risk_best,rev_best):
    a = ( risk - 0.01 ) * 100 / 0.08;
    b = ( rev  - 891  ) * 100 / 1455;

    f = (a*a + (100-b)*(100-b))^0.5;
else
    f = 145;
end

  out(1)=g1;
  out(2)=g2;
  out(3)=f;
  
end


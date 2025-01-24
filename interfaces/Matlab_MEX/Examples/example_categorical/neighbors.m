function [ Xout ] = neighbors ( X )  

t0 = X(1);
v0 = X(2);
t1 = X(3);
v1 = X(4);

t2 = ( 3 - t0 - t1 );

% neighbor #1:
Xout(1,:)=[ t2 ; v0 ; t1 ; v1 ];

% neighbor #2:
Xout(2,:)=[ t0 ; v0 ; t2 ; v1 ];

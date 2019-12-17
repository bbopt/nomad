close all
clear all
clc
disp('=========== EXPLORER =============');

%rand('twister',3)

! make -j 4
model = 'TYPE PRS DEGREE 2 RIDGE 0.01'
%model = 'TYPE RBF DISTANCE_TYPE NORM2_IS0 KERNEL_COEF 2 KERNEL_TYPE D2'
%model = 'TYPE KS DISTANCE_TYPE NORM2_IS0 KERNEL_COEF 10 KERNEL_TYPE D2'
%model = 'TYPE ENSEMBLE WEIGHT SELECT METRIC RMSECV'
model = 'TYPE ENSEMBLE WEIGHT SELECT METRIC OECV'
% model = 'TYPE KS'
% model = 'TYPE LOWESS'


%f = @(x) sum(x.^2,2);

N = 10;
P = 100;


M = 1;
X = randn(P,N);
Z = randn(P,M);

toBeDelData;
M = size(Z,2);

PP = 100;

dx = 0.1;
XX = zeros(PP,N);
for i=2:PP
    r = randn(1,N);
    r = dx*r/norm(r);
    XX(i,:) = XX(i-1,:)+r;
end


sgtelib_server_start(model,true)
sgtelib_server_ping;

sgtelib_server_newdata(X,Z);
[ZZmodel,std,ei,cdf] = sgtelib_server_predict(XX);
m = sgtelib_server_metric('RMSECV')

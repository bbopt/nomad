close all
clear all

% Start server
model = 'TYPE PRS DEGREE 2';
model = 'TYPE ENSEMBLE METRIC LINV WEIGHTS OPTIM'
sgtelib_server_start('TYPE PRS DEGREE 2',true)

% Test if server is ok and ready
sgtelib_server_ping;

% Build data points
N = 2;
M = 1;
P = 100;
X = rand(P,N);
Z = rand(P,M);

% Feed server
sgtelib_server_newdata(X,Z);

% Prediction points
PXX = 1000;
XX = rand(PXX,N);

%Prediction
[ZZ,std,ei,cdf] = sgtelib_server_predict(XX);

% Plot
figure; hold on;
plot3(X(:,1),X(:,2),Z,'*k');
plot3(XX(:,1),XX(:,2),ZZ,'or')

% Stop server
sgtelib_server_stop;
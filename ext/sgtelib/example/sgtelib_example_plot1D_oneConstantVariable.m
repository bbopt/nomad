close all
clear all
clc
disp('=========== EXPLORER =============');

! make -j 4

model = 'TYPE LOWESS DEGREE 2 KERNEL_SHAPE OPTIM METRIC OECV'

f = @(x) cos(4*x(:,1));



% Number of data points
P = 50;
% Data points
X = rand(P,2);
X(:,2) = 1;
Z = f(X);

Z(Z>prctile(Z,70)) = 1e20;




% Create model
PP = 50;
[x1,x2] = meshgrid( linspace(0,1,PP) , linspace(0,1,PP) );
XX = [x1(:) x2(:)];


sgtelib_server_start(model,true)
sgtelib_server_ping;
sgtelib_server_newdata(X,Z);
[ZZ,std,ei,cdf] = sgtelib_server_predict(XX);



figure('color','w'); hold on;
surf(x1,x2,reshape(ZZ,PP,PP));
plot3(X(:,1),X(:,2),Z,'ko','linewidth',3,'markersize',5);


set(gca,'view',[-36 66]);
zlim(bounds(ZZ));
xlabel('x');
ylabel('y');

sgtelib_server_stop;

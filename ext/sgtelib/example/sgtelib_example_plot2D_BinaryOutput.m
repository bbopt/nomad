close all
clear all
clc
disp('=========== EXPLORER =============');

! make -j 4

model = 'TYPE LOWESS DEGREE 1 KERNEL_SHAPE OPTIM METRIC OECV'

f = @(x) cos(sqrt(4*x(:,1))).*cos(3*x(:,2))+(x(:,1)>0.5)+(x(:,2)>0.5);



% Number of data points
P = 50;
% Data points
X = rand(P,2);
Z = f(X);
Z = rand + rand*(Z>median(Z));


% Create model
PP = 50;
[x1,x2] = meshgrid( linspace(0,1,PP) , linspace(0,1,PP) );
XX = [x1(:) x2(:)];


sgtelib_server_start(model,true)
sgtelib_server_ping;
sgtelib_server_newdata(X,Z);
[ZZ,std,ei,cdf] = sgtelib_server_predict(XX);



figure('color','w'); hold on;
imagesc(x1(1,:),x2(:,1),reshape(ZZ,PP,PP));
plot(X(:,1),X(:,2),'ko','linewidth',3,'markersize',5);

axis([0 1 0 1]);
xlabel('x');
ylabel('y');

sgtelib_server_stop;

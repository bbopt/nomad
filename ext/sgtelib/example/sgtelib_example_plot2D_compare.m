close all
clear all
clc
disp('=========== EXPLORER =============');

! make -j 4
model = 'TYPE LOWESS DISTANCE OPTIM KERNEL_SHAPE OPTIM KERNEL_TYPE OPTIM RIDGE OPTIM METRIC OECV'

model = 'TYPE KRIGING DISTANCE NORM2_IS0 METRIC OECV'

% creation of model
sgtelib_server_start(model,true)
sgtelib_server_ping;



% Test function
%f = @(x) cos(sqrt(0.5*x(:,1))).*cos(0.1*x(:,2));

f = @(x) 40*(x(:,1)==0)+...
         35*(x(:,2)==0)+...
         35*cos(x(:,1)/4)+...
         30*sin(x(:,2)/9);
  

% Dimension
N = 2;
% Number of sampling points along each dimension
PP1 = 25;
PP2 = 20;
% Bounds along each dimension
x1min = 0;
x1max = +30;
x2min = 0;
x2max = +30;

% Create data
scale_x1 = linspace(x1min,x1max,PP1)
scale_x2 = linspace(x2min,x2max,PP2)
[x1,x2] = meshgrid( scale_x1 , scale_x2 );
XX = [x1(:) x2(:)];
PP = (PP1)*(PP2);
ZZ = f(XX);


% Selection of some training points
P = 30;
i = randperm(PP);
i = i(1:P);
X = XX(i,:);
Z = ZZ(i,:);


% Sending data
sgtelib_server_newdata(X,Z);

% Prediction
[ZZ,std,ei,cdf] = sgtelib_server_predict(XX);



figure
subplot(2,2,1);hold on;
plot3(X(:,1),X(:,2),Z,'k.');
surf(x1,x2,reshape(f(XX),PP2,PP1));
set(gca,'view',[-37.5000   30.0000]);
set(gca,'ydir','normal')
xlabel('Real Z');
xlim([min(scale_x1) max(scale_x1)]);
ylim([min(scale_x2) max(scale_x2)]);
%zlim([-2 +2]);


subplot(2,2,2);hold on;
plot3(X(:,1),X(:,2),Z,'k.');
surf(x1,x2,reshape(ZZ,PP2,PP1));
set(gca,'view',[-37.5000   30.0000]);
xlabel('Zh');
xlim([min(scale_x1) max(scale_x1)]);
ylim([min(scale_x2) max(scale_x2)]);
%zlim([-2 +2]);

subplot(2,2,3);hold on;
imagesc(scale_x1,scale_x2,reshape(f(XX),PP2,PP1));
set(gca,'ydir','normal');
plot(X(:,1),X(:,2),'k.');
axis([min(scale_x1) max(scale_x1) min(scale_x2) max(scale_x2)])

subplot(2,2,4);hold on;
imagesc(scale_x1,scale_x2,reshape(ZZ,PP2,PP1));
set(gca,'ydir','normal');
plot(X(:,1),X(:,2),'k.');
axis([min(scale_x1) max(scale_x1) min(scale_x2) max(scale_x2)])



disp('metric:')
sgtelib_server_metric('AOECV')

sgtelib_server_stop;

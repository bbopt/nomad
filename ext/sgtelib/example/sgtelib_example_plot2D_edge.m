close all
clear all
clc
disp('=========== EXPLORER =============');

! make -j 4

%model = 'TYPE PRS'

model = 'TYPE ENSEMBLE DISTANCE OPTIM WEIGHT SELECT '

model = 'TYPE KRIGING DISTANCE OPTIM METRIC OECV'


%model = 'TYPE ENSEMBLE WEIGHT SELECT METRIC OECV DISTANCE OPTIM PRESET SUPER1'
%model = 'TYPE LOWESS KERNEL_TYPE OPTIM DISTANCE OPTIM'
%model = 'TYPE KS DISTANCE IS0'
%model = 'TYPE ENSEMBLE WEIGHT SELECT METRIC OECV DISTANCE IS0 PRESET SUPER1'


f = @(x) 40*(x(:,1)==0)+...
         35*(x(:,2)==0)+...
         35*cos(x(:,1)/4)+...
         30*sin(x(:,2)/9);
  

% Data points
x1max = 30;
x2max = 35;


p1 = 5;
p2 = 12;
p = 10;
X = [
     zeros(p1,1)              round(rand(p1,1)*x2max) ;...
     round(rand(p2,1 )*x1max) zeros(p2,1) ;...
     round(rand(p,1)*x1max)   round(rand(p,1)*x2max) ];
X = unique(X,'rows');
Z = f(X);


sgtelib_server_start(model,true)
sgtelib_server_ping;
sgtelib_server_newdata(X,Z);

% Prediction
[x1,x2] = meshgrid( (0:x1max) , (0:x2max) );
XX = [x1(:) x2(:)];
[ZZ,std,ei,cdf] = sgtelib_server_predict(XX);
ZZ = reshape(ZZ,x2max+1,x1max+1);


figure('color','w'); hold on;
surf(x1,x2,ZZ,'facealpha',0.5);
plot3(X(:,1),X(:,2),Z,'ko','linewidth',3,'markersize',5);
dv = max(Z)-min(Z);
zlim([min(Z)-dv/3,max(Z)+dv/3]);
set(gca,'view',[116 50]);
xlabel('x1');
ylabel('x2');

sgtelib_server_stop;

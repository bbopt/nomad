close all
%clear all
clc
disp('=========== EXPLORER =============');


model = 'TYPE ENSEMBLE WEIGHT SELECT METRIC OECV DISTANCE OPTIM'

model = 'TYPE KS DISTANCE OPTIM'
model = 'TYPE ENSEMBLE WEIGHT SELECT METRIC OECV DISTANCE OPTIM'

%model = 'TYPE PRS_CAT DEGREE 3'
%model = 'TYPE RBF DISTANCE OPTIM'

model = 'TYPE LOWESS DISTANCE OPTIM'


f = @(x)     x(:,2).^2   .* x(:,1) .* ( 2*mod(x(:,1),2)-1 )

% f = @(x) 40*(x(:,1)==0)+...
%          35*(x(:,2)==0)+...
%          35*cos(x(:,1)/4)+...
%          30*sin(x(:,2)/9);

% Data points
x1max = 5;
x2max = 100;

p = 50;
X = [round(rand(p,1)*x1max)   round(rand(p,1)*x2max) ];
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

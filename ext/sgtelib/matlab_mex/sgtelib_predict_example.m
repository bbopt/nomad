close all
clear all
clc

%model = 'TYPE ENSEMBLE   WEIGHT OPTIM   METRIC OECV';
model = 'TYPE LOWESS SHAPE_COEF OPTIM';
%model = 'TYPE RBF';

load('example15D.mat');

[ZZ ] = sgtelib_predict(LHS_X1,Y,LHS_X2,model);

plot([0.5 1 ],[0.5 1]);
hold on
plot(ZZ,YY,'Marker','o','LineStyle','none');

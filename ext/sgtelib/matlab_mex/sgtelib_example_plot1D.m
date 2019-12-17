close all
clear all
clc

% model = 'TYPE ENSEMBLE   WEIGHT OPTIM   METRIC OECV';
model = 'TYPE LOWESS SHAPE_COEF OPTIM';

P = 10;
X = sort(randn(P,1));

Z = X.^2-X-3;

PP = 5;
XX = sort([linspace(min(X)-1,max(X)+1,PP)' ; X]);

[ZZ ] = sgtelib_predict(X,Z,XX,model);

figure;
hold on;
plot(X,Z,'o');
plot(XX,ZZ,'-');

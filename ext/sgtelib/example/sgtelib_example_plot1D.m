close all
clear all
clc

model = 'TYPE ENSEMBLE   WEIGHT OPTIM   METRIC OECV';
model = 'TYPE LOWESS';
sgtelib_server_start(model,true)
sgtelib_server_ping;


P = 100;
X = sort(randn(P,1));


Z(:,1) = X.^2-X-3;
Z(:,2) = sign(X).*cos(X) + 0.1*randn(size(X));
Z(:,3) = mod(round(X),2)-0.5;
M = size(Z,2);
sgtelib_server_newdata(X,Z);



PP = 1000;
XX = sort([linspace(min(X)-1,max(X)+1,PP)' ; X]);



figure;
hold on;
for i=1:M
    plot(X,Z(:,i),'o','color',get_color(i,M));
end



[ZZ,std,ei,cdf] = sgtelib_server_predict(XX);
disp('metrics...')
sgtelib_server_metric('RMSE')
sgtelib_server_metric('RMSECV')
sgtelib_server_metric('OECV')

for i=1:M
    plot(XX,ZZ(:,i),'-','color',get_color(i,M),'linewidth',2);
    plot(XX,ZZ(:,i)+std(:,i),'--','color',get_color(i,M));
    plot(XX,ZZ(:,i)-std(:,i),'--','color',get_color(i,M));
end

% plot(xlim,[0 0],'--k');
% disp('Get CV values');
% [Zh,Sh,Zv,Sv] = sgtelib_server_cv;
% for i=1:M
%     plot(X,Zv(:,i),'o','color',get_color(i,M));
%     for j=1:P
%         plot(X(j)*[1 1],Zv(j,i)+Sv(j,i)*[-1 +1],'-','color',get_color(i,M));
%     end
% end


% disp('Plot EFI');
% figure;
% efi = ei(:,1).*prod(cdf(:,2:end),2);
% efi = max(1e-16,efi);
% plot(XX,efi,'g');
% set(gca,'yscale','log')



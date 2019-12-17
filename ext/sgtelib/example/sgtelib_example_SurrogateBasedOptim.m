close all
clear all
clc
disp('=========== EXPLORER =============');

addpath('../Matlab_Server');

WRITE = false;

! make -j 4
model = 'TYPE KRIGING'

sgtelib_server_start(model)



disp('===========   START  =============');
sgtelib_server_ping;


f = @(x) sin(5*x) + 1.5*x + x.^2;
xmin = -2;
xmax = +1.5;
zmin = -3;
zmax = +5;

X = [ -2 ; -1 ; 0.5 ; 1.5 ];
Z = f(X);

XX = linspace(xmin,xmax,1000)';

ZZtrue = f(XX);



figure('color','w'); hold on;
pointer = [];
leg_txt = cell(0);
k = 0;

pointer(1) = plot(XX,ZZtrue,'--b');
leg_txt{1} = 'True function f(x) (unknown)';
legend(pointer,leg_txt,'location','northwest'); legend('boxoff');
axis([xmin xmax zmin zmax]);
xlabel('x')
drawnow
if WRITE
    export_fig(['SurrogateBasedOptim_A' sprintf('%02d',k) '.pdf'],'-pdf'); k=k+1;
end

pointer(2) = plot(X,Z,'ok');
leg_txt{2} = 'Observations';
legend(pointer,leg_txt,'location','northwest'); legend('boxoff');
axis([xmin xmax zmin zmax]);
drawnow
if WRITE
    export_fig(['SurrogateBasedOptim_A' sprintf('%02d',k) '.pdf'],'-pdf'); k=k+1;
end
set(legend,'box','off')


sgtelib_server_newdata(X,Z);

for ii=1:8

    [ZZ,std,ei,cdf] = sgtelib_server_predict(XX);
    pointer(3) = plot(XX,ZZ,'r');
    if ii==1
        leg_txt{3} = 'Surrogate model';
    else
        leg_txt{3} = 'Updated surrogate model';
    end
    legend(pointer,leg_txt,'location','northwest'); legend('boxoff');
    axis([xmin xmax zmin zmax]);
    drawnow
    if WRITE
        export_fig(['SurrogateBasedOptim_A' sprintf('%02d',k) '.pdf'],'-pdf'); k=k+1;
    end

    [ZZmin,imin] = min(ZZ);
    XXmin = XX(imin);




    pointer(4) = plot(XXmin,ZZmin,'or');
    leg_txt{4} = ['Candidate: x=' num2str(XXmin,3) ];
    legend(pointer,leg_txt,'location','northwest'); legend('boxoff');
    axis([xmin xmax zmin zmax]);
    drawnow
    if WRITE
        export_fig(['SurrogateBasedOptim_A' sprintf('%02d',k) '.pdf'],'-pdf'); k=k+1;
    end

    Xnew = XXmin;
    Znew = f(Xnew);
    pointer(5) = plot(Xnew,Znew,'ob');
    leg_txt{5} = ['True value of candidate: f(' num2str(XXmin,3) ') = ' num2str(Znew,3)];
    legend(pointer,leg_txt,'location','northwest'); legend('boxoff');
    axis([xmin xmax zmin zmax]);
    drawnow
    if WRITE
        export_fig(['SurrogateBasedOptim_A' sprintf('%02d',k) '.pdf'],'-pdf'); k=k+1;
    end

    delete(pointer(3:5));
    pointer(3:5) = [];
    leg_txt(3:5) = [];
    plot(Xnew,Znew,'ko');
    sgtelib_server_newdata(Xnew,Znew);

end

sgtelib_server_stop;

close all
%clear all
clc
disp('=========== EXPLORER =============');

addpath('../Matlab_Server');

! make -j 4

model=cell(0,1);

 model{end+1} = 'TYPE LOWESS SHAPE_COEF -1 PRESET DEN '
 model{end+1} = 'TYPE LOWESS SHAPE_COEF -1 PRESET DGN '

 
 texOptions = {'fontname','times','fontsize',14};



disp('===========   START  =============');


% GOOD!!



X = randn(50,1);
Z = sin(2*pi*X)+0.1*randn(size(X));

Z = cos(pi*X).*sign(X);

%figure('color','w','position',[1466         527         455         437]); hold on;
figure('color','w','position',[198    95   933   639]); hold on;
pointer(1) = plot(X,Z,'ok','linewidth',2);


NN = 1000;
zoom = 1.0;
XX = sort([linspace(zoom*min(X),zoom*max(X),NN)' ; X ]);
%XX = 0.5;

metric = inf*ones(length(model),2);

for i=1:length(model)
    disp('Model : ');
    model{i}
    disp('Init');
    sgtelib_server_start(model{i},true)
    sgtelib_server_ping;
    disp('Data');
    sgtelib_server_newdata(X,Z);
    disp('Predict');
    [ZZ,std,ei,cdf] = sgtelib_server_predict(XX);
    metric(i,1) = sgtelib_server_metric('RMSECV');
    metric(i,2) = sgtelib_server_metric('OECV');
    marker = '';%get_marker(i,length(model));
    pointer(i+1) = plot(XX,ZZ,['-' marker],'color',get_color(i+1,length(model)+1));
    sgtelib_server_stop;   
end




legend_txt = {'Training points'};
for i=1:length(model)
    name = model{i};
    name = strrep(name,'TYPE ENSEMBLE WEIGHT SELECT METRIC','Select');
    name = strrep(name,'RMSECV','PRESS');
    name = strrep(name,'KERNEL_COEF ',' $\lambda$=');
    name = strrep(name,'TYPE',' ');
    name = strrep(name,'DEGREE','Degree');
    name = strrep(name,'RIDGE 0.001','Ridge 1e-3');
    name = strrep(name,'RIDGE 0.0',' ');
    name = strrep(name,'SHAPE_COEF',' ');
    name = strrep(name,'PRESET',' ');
    name = strrep(name,'Degree 2',' ');
    name = cleanSpaces(name);
    legend_txt{i+1} = name;
end
legend(pointer,legend_txt,'location','best','interpreter','latex');
%shift_legend_size([-0.15 0 +0.1 0]);
%axis([0 +1 -1.5 +1.5])

xlabel('x',texOptions{:});
ylabel('y',texOptions{:});

set(gca,texOptions{:});

 output_dir = './';
 output_file = [output_dir 'plot_function_2.pdf'];
export_fig([output_file],'-pdf');

for i=1:length(model)
    disp(legend_txt{i+1});
    s = '     ';
    s = [s num2str(metric(i,1),8)];
    s(20) = ' ';
    s = [s num2str(metric(i,2),8)];
    disp(s);
    
end
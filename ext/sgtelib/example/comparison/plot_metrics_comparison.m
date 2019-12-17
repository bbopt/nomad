function plot_metrics_comparison(input_file,output_file)
% output_dir = '.';
% input_file = 'output_hartman6.txt';
% output_dir = '.';
% output_file = [output_dir 'metric_comparison_hartman6.pdf'];
%close all



texOptions = {'fontname','times','fontsize',14};





m_order = [1 2 3 4];


model_list = cell(0);
metric_list = cell(0);
tab = [];

fid = fopen(input_file);
line = fgetl(fid);
while ischar(line)

    if isempty(line)
        line = fgetl(fid);
        continue;
    end
    line = cleanSpaces(line);

    w = getWords(line);
    if isempty(w{end})
        w(end)=[];
    end

    if strcmp(w{1},'metrics')
        w(1) = [];
        vs = zeros(1,length(w));
        for i = 1:length(w)
            vs(1,i) = str2num(w{i});
        end
        tab = [tab ; vs];
        % disregard data and delete last model
        if max(vs)<-1
            model_list(end) = [];
            tab(end,:) = [];
        end
    elseif strcmp(w{1},'list_metrics')
        w(1) = [];
        metric_list = w;
        disp(metric_list);
    elseif strcmp(w{1},'output')
        % nothing
    elseif strcmp(w{1},'TYPE')
        disp('=================')
        line
       
        line = cleanSpaces(line);
        line
        line = strrep(line,'TYPE','');
        line = strrep(line,'DEGREE','');
        line = strrep(line,'KERNEL_COEF -1','OPT');
        line = strrep(line,'KERNEL_COEF','');
        line = strrep(line,'KERNEL_','');
        line = strrep(line,'D1','');
        line = strrep(line,'WEIGHT','');
        line = strrep(line,'METRIC','');
        line = strrep(line,'ENSEMBLE','Ens.');
        line = strrep(line,'OPTIM','Optim.');
        line = strrep(line,'SELECT','Select');
        line = strrep(line,'RIDGE 0.001','Ridge 1e-3');
        line = strrep(line,'RIDGE 0','');
        line = cleanSpaces(line);
        line = strrep(line,'RBFI I','RBFI PolyH. ');
        line = strrep(line,'PRESET',' ');
        line = cleanSpaces(line);
        model_list{end+1} = line;
    else
        disp(line)
        error('line not recognized');
    end

    line = fgetl(fid);
end

fclose(fid);



% Count
NM = length(metric_list);
NS = length(model_list);

% Reorder
metric_list = metric_list(m_order);
tab = tab(:,m_order);

% Build fancy names
metric_fancy_names = metric_list;
metric_fancy_names = strrep(metric_fancy_names,'LOO','RMSE');
metric_fancy_names = strrep(metric_fancy_names,'RMSECV','PRESS');
metric_fancy_names = strrep(metric_fancy_names,'OETP','OE');


xmin = zeros(1,NM);
xmax = zeros(1,NM);
for nm=1:NM
    e = tab(:,nm);
    xmin(nm) = min(e(e>0));
    if any(e==0)
        xmin(nm) = xmin(nm)/10;
    end
    xmax(nm) = max(e);
end

xmin = 10.^floor(log10(xmin));
xmax = 10.^ceil (log10(xmax));

% 
% 
% xmin(:) = 0.0;
% xmax(:) = 1;





MODEL_TYPES = {'PRS','KS','RBFI','LOWESS','Ens.'};

%MODEL_TYPES = {'PRS','KS','RBFI','RBFI PolyH'};


%figure('color','w','position',[959    19   932   945]);
%figure('color','w','position',[949   483   932   481]);
%figure('color','w','position',[793    14   648   770]);
%figure('color','w','position',[509   574 932   210]);
figure('color','w','position',[586   408   847   374]);
set(gcf,'units','normalized');

coladd = 2.5
height = 0.85;
width = 0.7;
for nm=1:NM
    
    % ORder
    t = tab(:,nm)
    [t,t] = sort(t);
    t(t) = (1:length(t))
    
    pos = [ (nm+coladd-1)/(NM+coladd) (1-height)/2 width/(NM+coladd) height];
    subplot('position',pos); hold on;
    for ns=1:NS
        w = getWords(model_list{ns});
        c = -1;
        for i=1:length(MODEL_TYPES)
            if ~isempty(strfind(model_list{ns},MODEL_TYPES{i}))
                c = i;
            end
        end
        if c==-1
            w{1}
            MODEL_TYPES
            error('not recognised...');
        end
        c = get_color(c+1,length(MODEL_TYPES)+1);
        c = 1-0.7*(1-c);
        fill([xmin(nm)*[1 1] tab(ns,nm)*[1 1]],ns+[-0.4 +0.4 +0.4 -0.4],c);
        text(0.75,ns,num2str(t(ns)));
        %text(0.75,ns, num2str(t(ns)));
    end


    set(gca,'units','normalized');
    title(metric_fancy_names{nm},texOptions{:});
    if nm==1
        disp('Modify ytick');
        set(gca,'ytick',(1:NS),'yticklabel',model_list,texOptions{:});
    else
        set(gca,'ytick',[],'yticklabel',[],texOptions{:});
    end
    %set(gca,'xscale','log');
    %xlim([xmin(nm) xmax(nm)]);
    %set(gca,'xtick',[xmin(nm) xmax(nm)]);
    set(gca,'ydir','reverse','xgrid','on','xminorgrid','off','box','on');
    ylim([0.5 NS+0.5]);
    
    
    

end

export_fig([output_file],'-pdf');













output_file = strrep(output_file,'.pdf','_correlation.tex');
fid = fopen(output_file,'w');
newline = char(10);

s = '\begin{tabular}{| c ||';
for i=1:NM
    s = [s ' c |'];
end
fwrite(fid,[s '}' newline]);
fwrite(fid,['\hline' newline]);
s = '';
for i=1:NM
    s = [s ' & ' metric_fancy_names{i} ];
end
fwrite(fid,[s ' \\' newline]);
fwrite(fid,['\hline\hline' newline]);
for i=1:NM
    mi = tab(:,i);
    fwrite(fid,[metric_fancy_names{i}  newline]);
    for j=1:NM
        if j>=i
            mj = tab(:,j);        
            c = corr([mi mj]);
            c = c(2,1);
            c = 0.01*round(c*100);
            s = num2str(c);
            if abs(c)~=1 && c~=0 && 10*c==round(10*c)
                s = [s '0'];
                if isempty(find(s=='.'))
                    error('PROBLEM...');
                end
            end
        else
            s = ' ';
        end
        if j==NM
            s = [s ' \\'];
        end
        fwrite(fid,[ ' & ' s newline]);
    end
    fwrite(fid,['\hline' newline]);
end
fwrite(fid,'\end{tabular}');
fclose(fid);
%edit(output_file);








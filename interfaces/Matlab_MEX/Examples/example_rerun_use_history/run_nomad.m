%% TEST IS NOT WORKING
%% TODO fix it if needed



global historyData useHistoryData

clc
x0 = [1;1];

useHistoryData =true;
historyFileInit = 'history.0.txt';

historyData=[];
if useHistoryData
    historyData = load(historyFileInit);
end

params = struct('display_degree','2','bb_output_type','OBJ EB','history_file','history.txt','max_bb_eval','10');
[x,fval] = nomadOpt(@bb,x0,[-10;-10],[10;10],params);

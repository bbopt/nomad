function eval = bb(x)

global historyData useHistoryData

dataAvailable = false;
if useHistoryData
   [C,index]=intersect(historyData(:,1:2),x','rows','stable');
   
   if size(index)~=0
       dataAvailable = true;
       
   end
end
   
if dataAvailable
    eval = historyData(index,3:4); 
else
    eval=[10*(x(2)-x(1)^2), 1 - x(1)];
end

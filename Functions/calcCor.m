function numCor = calcCor(data)

% data = sort(data(:,[2 3]),1);
data = data(:,[2 3]);
c = unique(data(:,1));
numCor = cat(2,c,zeros(length(c),1));
for cnd = 1:length(c) 
    curCnd = data(find(data == c(cnd)),:);
    sumCor = sum(curCnd(:,2));
    numCor(cnd,2) = sumCor;
end
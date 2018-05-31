function percentCor = calcCor(data)
% takes in correctness counts per crowding condition and returns accuracy
% percentages
% INPUT: nTrial x [crowding str, angle, correctness]

data = data(:,[2 4]);

% sum number of each condition
uniqueCnd = unique(data(:,1));

% pre-allocate array for input of correct answer counts
percentCor = cat(2,uniqueCnd,zeros(length(uniqueCnd),1));

for i = 1:length(uniqueCnd)
    curIdx = find(data(:,1) == uniqueCnd(i));
    numPerCnd = length(curIdx);
    curCnd = data(curIdx,:);
    percentCor(i,2) = sum(curCnd(:,2))/numPerCnd;
end
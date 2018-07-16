%% DATASET INFO
% 1 Bethany
% 2 Jianfei - 2x threshold
% 3 Dan - 2x threshold, mask between S1 and S2 instead of blank
% 4 Bethany - 2x threshold, mask between S1 and S2 instead of blank
% 5 Sydney - 2x threshold, mask between S1 and S2 instead of blank

%%
clear all
close all

% parentDir = '~/Bethany/paclab';
parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));
dataDir = sprintf('%s/Subject_folders',parentDir);
subj = sort(strsplit(ls(dataDir)));
subj = subj(6:end);

for s = 1:length(subj)
    curSubj = subj{s};
    subjNum = str2num(curSubj(2:3));
    subjDir = sprintf('%s/%s',dataDir,curSubj);
    
    dataThresh = load(sprintf('%s/%dblock1_threshold_all.mat',subjDir,subjNum));
    thresholds(s,:) = dataThresh.thresholds;
    
    blockFiles = mySubFiles(subjDir,'noAction_main_all.mat',9);
    for b = 1:length(blockFiles)
        data = load(sprintf('%s/%s',subjDir,blockFiles{b}));
        targID = data.targID;
        rspKey = data.rspKey;
        cndList = data.cndList;
        catchIdx = data.catchIdx;
        all = 1:length(targID);
        nonCatchIdx = setdiff(all,catchIdx);
        actualTargID(1+(b-1)*length(nonCatchIdx):b*length(nonCatchIdx)) = targID(nonCatchIdx);
        actualRspKey(1+(b-1)*length(nonCatchIdx):b*length(nonCatchIdx)) = rspKey(nonCatchIdx);
    end
    
    hit = 0; miss = 0; falseAlarm = 0; corReject = 0;
    for i = 1:length(actualTargID)
        if actualTargID(i) && actualRspKey(i)
            hit = hit + 1;
        elseif actualTargID(i) && ~actualRspKey(i)
            miss = miss + 1;
        elseif ~actualTargID(i) && actualRspKey(i)
            falseAlarm = falseAlarm + 1;
        elseif ~actualTargID(i) && ~actualRspKey(i)
            corReject = corReject + 1;
        end
    end
    
    yes(s) = sum(actualTargID); % shoudl always be 144 / need to fix catch trial randomization
    no(s) = length(actualTargID) - sum(actualTargID); 
    H(s) = hit/yes(s);
    F(s) = falseAlarm/no(s);
    dp(s) = norminv(H(s)) - norminv(F(s));
    c = -0.5*(norminv(H)+ norminv(F));
end
%%

gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};
figure
bar(thresholds');
set(gca,gcaOpts{:})
ylim([0 8])

%%
grandMean = mean(stairmean,1); 
cndNames = {'Subj 1' 'Subj 2' 'Subj 3' 'Subj 4'};
gcaOpts = {'XTick',1:4,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

% REORGANIZING & COLLAPSING DATA
sepCrwd = [mean(finalRev(:,1:2),2), mean(finalRev(:,3:4),2)];
for i = 1:size(final6Rev,3)
    curNC = final6Rev(:,1:2,i);
    curC = final6Rev(:,3:4,i);
    avgData(i,:) = [mean2(curNC) mean2(curC)];
    semData(i,:) = [std2(curNC) std2(curC)]/sqrt(length(curNC)*2);
end

%% PLOT THRESHOLDS
close all
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
h = bar(avgData'); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
h(4).FaceColor = [28 15 142]/255;
if dataset == 3
    hold on;
    ngroups = size(avgData, 2);
    nbars = size(avgData, 1);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    for i = 1:nbars
        x = (1:ngroups) - (groupwidth/2) + (2*i-1) * (groupwidth / (2*nbars));
        errorbar(x, avgData(i,:), semData(i,:), 'k', 'linestyle', 'none');
    end
    legend({subj{:}, 'SEM'},'AutoUpdate','off','location','northeast');
else
    legend({subj{:}},'AutoUpdate','off','location','northeast');
end
set(gca,gcaOpts{:})
title('Avg of last 6 reversal values, avged over up/down staircases')
xlabel('Crowding condition')
ylabel('Avg threshold (deg)')
ylim([0 20])


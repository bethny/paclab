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

noCrwdTargID = []; noCrwdRspKey = []; crwdTargID = []; crwdRspKey = [];

for s = 1:length(subj)
    curSubj = subj{s};
    subjNum = str2num(curSubj(2:3));
    subjDir = sprintf('%s/%s',dataDir,curSubj);
    blockFiles = mySubFiles(subjDir,'noAction_main_all.mat',9);
    
    dataThresh = load(sprintf('%s/%dblock1_threshold_all.mat',subjDir,subjNum));
    thresholds(s,:) = dataThresh.thresholds;
    
    for b = 1:length(blockFiles)
        data = load(sprintf('%s/%s',subjDir,blockFiles{b}));
        targID = data.targID;
        rspKey = data.rspKey;
        cndList = data.cndList;
        catchIdx = data.catchIdx;
        all = 1:length(targID);
        nonCatchIdx = setdiff(all,catchIdx); % SHOULD BE 144 TRIALS ALWAYS 
        
        actualTargID = targID(nonCatchIdx);
        actualRspKey = rspKey(nonCatchIdx);
        cndList = cndList(nonCatchIdx);
        noCrwd_idx = find(cndList == 5 | cndList == 6 | cndList == 7 | cndList == 8); % should be 48 pre-correction
        crwd_idx = find(cndList == 1 | cndList == 2 | cndList == 3 | cndList == 4 ...
             | cndList == 9 | cndList == 10 | cndList == 11 | cndList == 12); % should be 96 pre-correction
        
        x = actualTargID(noCrwd_idx);
        y = actualRspKey(noCrwd_idx);
        z = actualTargID(crwd_idx);
        a = actualRspKey(crwd_idx);
        
        noCrwdTargID = cat(2,noCrwdTargID,x);
        noCrwdRspKey = cat(2,noCrwdRspKey,y);
        crwdTargID = cat(2,crwdTargID,z);
        crwdRspKey = cat(2,crwdRspKey,a);
    end
    
    % FOR NO CROWDING CONDITIONS
    NC_hit = 0; NC_miss = 0; NC_falseAlarm = 0; NC_corReject = 0;
    for i = 1:length(noCrwdTargID)
        if noCrwdTargID(i) && noCrwdRspKey(i)
            NC_hit = NC_hit + 1;
        elseif noCrwdTargID(i) && ~noCrwdRspKey(i)
            NC_miss = NC_miss + 1;
        elseif ~noCrwdTargID(i) && noCrwdRspKey(i)
            NC_falseAlarm = NC_falseAlarm + 1;
        elseif ~noCrwdTargID(i) && ~noCrwdRspKey(i)
            NC_corReject = NC_corReject + 1;
        end
    end
    
    NC_yes(s) = sum(noCrwdTargID); % shoudl always be 144 / need to fix catch trial randomization
    NC_no(s) = length(noCrwdTargID) - sum(noCrwdTargID); 
    NC_H(s) = NC_hit/NC_yes(s);
    NC_F(s) = NC_falseAlarm/NC_no(s);
    NC_dp(s) = norminv(NC_H(s)) - norminv(NC_F(s));
    NC_c = -0.5*(norminv(NC_H)+ norminv(NC_F));
    
    % FOR CROWDING CONDITIONS
    C_hit = 0; C_miss = 0; C_falseAlarm = 0; C_corReject = 0;
    for i = 1:length(crwdTargID)
        if crwdTargID(i) && crwdRspKey(i)
            C_hit = C_hit + 1;
        elseif crwdTargID(i) && ~crwdRspKey(i)
            C_miss = C_miss + 1;
        elseif ~crwdTargID(i) && crwdRspKey(i)
            C_falseAlarm = C_falseAlarm + 1;
        elseif ~crwdTargID(i) && ~crwdRspKey(i)
            C_corReject = C_corReject + 1;
        end
    end
    
    C_yes(s) = sum(crwdTargID); % shoudl always be 144 / need to fix catch trial randomization
    C_no(s) = length(crwdTargID) - sum(crwdTargID); 
    C_H(s) = C_hit/C_yes(s);
    C_F(s) = C_falseAlarm/C_no(s);
    C_dp(s) = norminv(C_H(s)) - norminv(C_F(s));
    C_c = -0.5*(norminv(C_H)+ norminv(C_F));
    
end

%% PLOT THRESHOLDS
close all
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
x1 = thresholds(:,1)';
x1(length(thresholds) + 1) = mean(thresholds(:,1));
x2 = thresholds(:,2)';
x2(length(thresholds) + 1) = mean(thresholds(:,2));

h = bar([x1;x2]);
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
h(4).FaceColor = [28 15 142]/255;
h(5).FaceColor = [255 255 255]/255;
h(6).FaceColor = [0 0 0]/255;
legend({subj{:},'average'},'AutoUpdate','off','location','northwest');
set(gca,gcaOpts{:})
title('Avg of last 6 reversal values, avged over up/down staircases')
xlabel('Crowding condition')
ylabel('Avg threshold (deg)')
ylim([0 10])

%% PLOT D PRIMES
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
h = bar([NC_dp, mean(NC_dp); C_dp, mean(C_dp)]); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
h(4).FaceColor = [28 15 142]/255;
h(5).FaceColor = [255 255 255]/255;
h(6).FaceColor = [0 0 0]/255;
legend({subj{:}, 'average'},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('D-prime values')
xlabel('Crowding condition')
ylabel('D-prime')
ylim([0 3])

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
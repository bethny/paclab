%% DATASET INFO
% 1 Bethany
% 2 Jianfei - 2x threshold
% 3 Dan - 2x threshold, mask between S1 and S2 instead of blank
% 4 Bethany - 2x threshold, mask between S1 and S2 instead of blank
% 5 Sydney - 2x threshold, mask between S1 and S2 instead of blank

%%
clear all
clear all
close all

% WHICH DATASET?
dataset = 1; % 0 = main exp no grasp, 1 = main exp with grasp

% parentDir = '~/Bethany/paclab';
parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));
dataDir = sprintf('%s/Subject_folders',parentDir);
if ~dataset
    dataDir = sprintf('%s/Pilot_NoGrasp',dataDir);
end
subj = sort(strsplit(ls(dataDir)));
subj = subj(7:end);

NC_TargID = []; NC_RspKey = []; C_TargID = []; C_RspKey = [];
NC_TargID_key = []; NC_RspKey_key = []; C_TargID_key = []; C_RspKey_key = [];

for s = 1:length(subj)
    curSubj = subj{s};
    subjNum = str2num(curSubj(2:3));
    subjDir = sprintf('%s/%s',dataDir,curSubj);
    blockFiles = mySubFiles(subjDir,'all',13);
    
    dataThresh = load(blockFiles{1});
    thresholds(s,:) = dataThresh.thresholds;
    
    for b = 2:length(blockFiles)
        data = load(blockFiles{b});
        targID = data.targID;
        rspKey = data.rspKey;
        cndList = data.cndList;
        catchIdx = data.catchIdx;
        all = 1:length(targID);
        nonCatchIdx = setdiff(all,catchIdx); % LENGTH = 96 TRIALS 
        
        actualTargID = targID(nonCatchIdx);
        actualRspKey = rspKey(nonCatchIdx);
        cndList = cndList(nonCatchIdx);
        NC_idx = find(cndList == 5 | cndList == 6 | cndList == 7 | cndList == 8); % length = 48 / non crowd
        C_idx = find(cndList == 1 | cndList == 2 | cndList == 3 | cndList == 4 ...
             | cndList == 9 | cndList == 10 | cndList == 11 | cndList == 12); % length = 48 / crowd
        
        x = actualTargID(NC_idx);
        y = actualRspKey(NC_idx);
        z = actualTargID(C_idx);
        a = actualRspKey(C_idx);
        
        if ~data.grasping
            NC_TargID_key = cat(2,NC_TargID_key,x);
            NC_RspKey_key = cat(2,NC_RspKey_key,y);
            C_TargID_key = cat(2,C_TargID_key,z);
            C_RspKey_key = cat(2,C_RspKey_key,a);
        else
            NC_TargID = cat(2,NC_TargID,x);
            NC_RspKey = cat(2,NC_RspKey,y);
            C_TargID = cat(2,C_TargID,z);
            C_RspKey = cat(2,C_RspKey,a);
        end
    end
    
    yes = sum(NC_TargID); % should always be 48
    no = length(NC_TargID) - sum(NC_TargID); 
    
    % FOR NO CROWDING CONDITIONS
    % KEYPRESS ONLY
    NC_hit_key = 0; NC_miss_key = 0; NC_falseAlarm_key = 0; NC_corReject_key = 0;
    for i = 1:length(NC_TargID_key)
        if NC_TargID_key(i) && NC_RspKey_key(i)
            NC_hit_key = NC_hit_key + 1;
        elseif NC_TargID_key(i) && ~NC_RspKey_key(i)
            NC_miss_key = NC_miss_key + 1;
        elseif ~NC_TargID_key(i) && NC_RspKey_key(i)
            NC_falseAlarm_key = NC_falseAlarm_key + 1;
        elseif ~NC_TargID_key(i) && ~NC_RspKey_key(i)
            NC_corReject_key = NC_corReject_key + 1;
        end
    end
    
    % GRASP
    NC_hit = 0; NC_miss = 0; NC_falseAlarm = 0; NC_corReject = 0;
    for i = 1:length(NC_TargID)
        if NC_TargID(i) && NC_RspKey(i)
            NC_hit = NC_hit + 1;
        elseif NC_TargID(i) && ~NC_RspKey(i)
            NC_miss = NC_miss + 1;
        elseif ~NC_TargID(i) && NC_RspKey(i)
            NC_falseAlarm = NC_falseAlarm + 1;
        elseif ~NC_TargID(i) && ~NC_RspKey(i)
            NC_corReject = NC_corReject + 1;
        end
    end    
       
    NC_H_key(s) = NC_hit_key/yes;
    NC_F_key(s) = NC_falseAlarm_key/no;
    if ~NC_F_key(s)
        NC_F_key(s) = 0.01;
    end
    NC_dp_key(s) = norminv(NC_H_key(s)) - norminv(NC_F_key(s));
    
    NC_H(s) = NC_hit/yes;
    NC_F(s) = NC_falseAlarm/no;
    if ~NC_F(s)
        NC_F(s) = 0.01;
    end
    NC_dp(s) = norminv(NC_H(s)) - norminv(NC_F(s));
  
    % FOR CROWDING CONDITIONS
    % KEYPRESS ONLY
    C_hit_key = 0; C_miss_key = 0; C_falseAlarm_key = 0; C_corReject_key = 0;
    for i = 1:length(C_TargID_key)
        if C_TargID_key(i) && C_RspKey_key(i)
            C_hit_key = C_hit_key + 1;
        elseif C_TargID_key(i) && ~C_RspKey_key(i)
            C_miss_key = C_miss_key + 1;
        elseif ~C_TargID_key(i) && C_RspKey_key(i)
            C_falseAlarm_key = C_falseAlarm_key + 1;
        elseif ~C_TargID_key(i) && ~C_RspKey_key(i)
            C_corReject_key = C_corReject_key + 1;
        end
    end
    
    % GRASP
    C_hit = 0; C_miss = 0; C_falseAlarm = 0; C_corReject = 0;
    for i = 1:length(C_TargID)
        if C_TargID(i) && C_RspKey(i)
            C_hit = C_hit + 1;
        elseif C_TargID(i) && ~C_RspKey(i)
            C_miss = C_miss + 1;
        elseif ~C_TargID(i) && C_RspKey(i)
            C_falseAlarm = C_falseAlarm + 1;
        elseif ~C_TargID(i) && ~C_RspKey(i)
            C_corReject = C_corReject + 1;
        end
    end   
    
    C_H_key(s) = C_hit_key/yes;
    C_F_key(s) = C_falseAlarm_key/no;
    C_dp_key(s) = norminv(C_H_key(s)) - norminv(C_F_key(s));
    
    C_H(s) = C_hit/yes;
    C_F(s) = C_falseAlarm/no;
    C_dp(s) = norminv(C_H(s)) - norminv(C_F(s));
    
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

%% PLOT D PRIMES / WITH GRASP
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
h = bar([NC_dp_key; NC_dp; C_dp_key; C_dp]); 
h(1).FaceColor = [135,205,215]/255;
% h(2).FaceColor = [240,59,37]/255;
% h(3).FaceColor = [250,185,219]/255;
% h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({'BH'},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('D-prime values')
xlabel('Condition')
ylabel('D-prime')
ylim([0 4])


%% PLOT D PRIMES / NO GRASP
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
h = bar([NC_dp_key, mean(NC_dp_key); C_dp_key, mean(C_dp_key)]); 
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
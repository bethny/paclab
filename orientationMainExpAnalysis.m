%% DATASET INFO
% MAIN EXP WITHOUT GRASP, FIRST ATTEMPT
% 1 Bethany
% 2 Jianfei - 2x threshold
% 3 Dan - 2x threshold, mask between S1 and S2 instead of blank
% 4 Bethany - 2x threshold, mask between S1 and S2 instead of blank
% 5 Sydney - 2x threshold, mask between S1 and S2 instead of blank

% MAIN EXP, FULL VERSION
% 1 Bethany - weird thresholds, 1.5x multiplier 
% 2 Dan - 1.25x multiplier
% 3 Will - 1.25x multiplier

% calculate reaction time after s2 onset + initiation latency for movement

%%
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
nSubj = length(subj);

% PREALLOCATION
NC_hit_key = zeros(1,nSubj); NC_miss_key = zeros(1,nSubj); NC_falseAlarm_key = zeros(1,nSubj); NC_corReject_key = zeros(1,nSubj);
NC_H_key = zeros(1,nSubj); NC_F_key = zeros(1,nSubj); NC_dp_key = zeros(1,nSubj);
C_hit_key = zeros(1,nSubj); C_miss_key = zeros(1,nSubj); C_falseAlarm_key = zeros(1,nSubj); C_corReject_key = zeros(1,nSubj);
C_H_key = zeros(1,nSubj); C_F_key = zeros(1,nSubj); C_dp_key = zeros(1,nSubj);
NC_hit = zeros(1,nSubj); NC_miss = zeros(1,nSubj); NC_falseAlarm = zeros(1,nSubj); NC_corReject = zeros(1,nSubj);
NC_H = zeros(1,nSubj); NC_F = zeros(1,nSubj); NC_dp = zeros(1,nSubj);
C_hit = zeros(1,nSubj); C_miss = zeros(1,nSubj); C_falseAlarm = zeros(1,nSubj); C_corReject = zeros(1,nSubj);
C_H = zeros(1,nSubj); C_F = zeros(1,nSubj); C_dp = zeros(1,nSubj);
nInvalidKey = zeros(nSubj,2); percAccGrasp = zeros(nSubj,2);
diff = zeros(1,nSubj); same = zeros(1,nSubj); catch_diff = zeros(1,nSubj); catch_same = zeros(1,nSubj); 

for s = 1:nSubj
    NC_TargID = []; NC_RspKey = []; C_TargID = []; C_RspKey = [];
    NC_TargID_key = []; NC_RspKey_key = []; C_TargID_key = []; C_RspKey_key = [];
    catch_NC_TargID = []; catch_NC_RspKey = []; catch_C_TargID = []; catch_C_RspKey = [];
    catch_NC_TargID_key = []; catch_NC_RspKey_key = []; catch_C_TargID_key = []; catch_C_RspKey_key = [];
    
    curSubj = subj{s};
    subjNum = str2num(curSubj(2:3));
    subjDir = sprintf('%s/%s',dataDir,curSubj);
    blockFiles = mySubFiles(subjDir,'all',13);
    
    dataThresh = load(blockFiles{1});
    thresholds(s,:) = dataThresh.thresholds;
    
    idx = 1; grIdx = 1;
    for b = 2:length(blockFiles)
        data = load(blockFiles{b});
        targID = data.targID;
        rspKey = data.rspKey;
        allCndList = data.cndList;
        catchIdx = data.catchIdx;
        nonCatchIdx = data.nonCatchIdx;
        nActual = data.trialNumber;
        nCatch = data.nCatch;
        
        acc = data.acc;
        S1_onset_time = data.S1_onset_time;
        S2_onset_time = data.S2_onset_time;
        
        rt = data.rt;
        
        % RT ANALYSIS (for all trials)
        
        if ~data.grasping
           nInvalidKey(s,idx) = length(find(rt<(.1 + .3 + .1 + .1))); 
           idx = idx + 1;
        else
           SOT_data = data.SOT_data;
           percAccGrasp(s,grIdx) = sum(data.acc_grasp)/length(data.acc_grasp);
%            initLatency = SOT_data(end) - S1_onset_time 
            grIdx = grIdx + 1;
        end
        
        allTrials = 1:length(targID); % 128 trials, per block
        
        actualTargID = targID(nonCatchIdx); % len = 96
        actualRspKey = rspKey(nonCatchIdx); % len = 96
        cndList = allCndList(nonCatchIdx); % len = 96
        
        catchTargID = targID(catchIdx); % len = 32
        catchRspKey = rspKey(catchIdx); % len = 32
        catchCndList = allCndList(catchIdx); % len = 32 
        
        NC_idx = find(cndList == 5 | cndList == 6 | cndList == 7 | cndList == 8); % length = 48 / non crowd
        C_idx = find(cndList == 1 | cndList == 2 | cndList == 3 | cndList == 4 ...
             | cndList == 9 | cndList == 10 | cndList == 11 | cndList == 12); % length = 48 / crowd
        catch_NC_idx = find(catchCndList == 5 | catchCndList == 6 | catchCndList == 7 | catchCndList == 8); % length = 16 / NC
        catch_C_idx = find(catchCndList == 1 | catchCndList == 2 | catchCndList == 3 | catchCndList == 4 ...
             | catchCndList == 9 | catchCndList == 10 | catchCndList == 11 | catchCndList == 12); % length = 16 / crowd
        
        x = actualTargID(NC_idx); % NONCROWDING TARGET ID; length = 48 (per block)
        y = actualRspKey(NC_idx); % NONCROWDING RESP KEY; length = 48 (per block)
        z = actualTargID(C_idx); % CROWDING TARGET ID; length = 48 (per block)
        a = actualRspKey(C_idx); % CROWDING RESP KEY; length = 48 (per block)
        
        r = catchTargID(catch_NC_idx); % CATCH NONCROWDING TARGET ID; length = 16 (per block)
        t = catchRspKey(catch_NC_idx); % CATCH NONCROWDING RESP KEY; length = 16 (per block)
        u = catchTargID(catch_C_idx); % CATCH CROWDING TARGET ID; length = 16 (per block)
        v = catchRspKey(catch_C_idx); % CATCH CROWDING RESP KEY; length = 16 (per block)
        
        if ~data.grasping
            NC_TargID_key = cat(2,NC_TargID_key,x); % len = 96, because 48 per block * 2 blocks
            NC_RspKey_key = cat(2,NC_RspKey_key,y); % len = 96, because 48 per block * 2 blocks
            C_TargID_key = cat(2,C_TargID_key,z); % len = 96, because 48 per block * 2 blocks
            C_RspKey_key = cat(2,C_RspKey_key,a); % len = 96, because 48 per block * 2 blocks
            
            catch_NC_TargID_key = cat(2,catch_NC_TargID_key,r); % len = 32, because 16 per block * 2 blocks
            catch_NC_RspKey_key = cat(2,catch_NC_RspKey_key,t); % len = 32, because 16 per block * 2 blocks
            catch_C_TargID_key = cat(2,catch_C_TargID_key,u); % len = 32, because 16 per block * 2 blocks
            catch_C_RspKey_key = cat(2,catch_C_RspKey_key,v); % len = 32, because 16 per block * 2 blocks
        else
            NC_TargID = cat(2,NC_TargID,x); % len = 96, because 48 per block * 2 blocks
            NC_RspKey = cat(2,NC_RspKey,y); % len = 96, because 48 per block * 2 blocks
            C_TargID = cat(2,C_TargID,z); % len = 96, because 48 per block * 2 blocks
            C_RspKey = cat(2,C_RspKey,a); % len = 96, because 48 per block * 2 blocks
            
            catch_NC_TargID = cat(2,catch_NC_TargID,r); % len = 32, because 16 per block * 2 blocks
            catch_NC_RspKey = cat(2,catch_NC_RspKey,t); % len = 32, because 16 per block * 2 blocks
            catch_C_TargID = cat(2,catch_C_TargID,u); % len = 32, because 16 per block * 2 blocks
            catch_C_RspKey = cat(2,catch_C_RspKey,v); % len = 32, because 16 per block * 2 blocks
        end
    end
    
    diff(s) = sum(NC_TargID); % = 48, because 48 total DIFF trials per task type (grasp/key) per subj; 24 * 2
    same(s) = length(NC_TargID) - sum(NC_TargID); % = 48, same reasoning ^
    
    catch_diff(s) = sum(catch_NC_TargID);% = 16, because 16 total DIFF trials per task type per subj; 8 * 2
    catch_same(s) = length(catch_NC_TargID) - sum(catch_NC_TargID); % = 16, same reasoning ^
    
    % FOR NO CROWDING CONDITIONS %%%%%%%%%
    % KEYPRESS ONLY %%%%
    NC_hit_key(s) = 0; NC_miss_key(s) = 0; NC_falseAlarm_key(s) = 0; NC_corReject_key(s) = 0;
    for i = 1:length(NC_TargID_key)
        if NC_TargID_key(i) && NC_RspKey_key(i)
            NC_hit_key(s) = NC_hit_key(s) + 1;
        elseif NC_TargID_key(i) && ~NC_RspKey_key(i)
            NC_miss_key(s) = NC_miss_key(s) + 1;
        elseif ~NC_TargID_key(i) && NC_RspKey_key(i)
            NC_falseAlarm_key(s) = NC_falseAlarm_key(s) + 1;
        elseif ~NC_TargID_key(i) && ~NC_RspKey_key(i)
            NC_corReject_key(s) = NC_corReject_key(s) + 1;
        end
    end
    
    NC_H_key(s) = NC_hit_key(s)/diff(s);
    if ~NC_H_key(s)
        NC_H_key(s) = 0.01;
    end
    NC_F_key(s) = NC_falseAlarm_key(s)/same(s);
    if ~NC_F_key(s)
        NC_F_key(s) = 0.01;
    end
    NC_dp_key(s) = norminv(NC_H_key(s)) - norminv(NC_F_key(s));
    
    % GRASP %%%%
    NC_hit(s) = 0; NC_miss(s) = 0; NC_falseAlarm(s) = 0; NC_corReject(s) = 0;
    for i = 1:length(NC_TargID)
        if NC_TargID(i) && NC_RspKey(i)
            NC_hit(s) = NC_hit(s) + 1;
        elseif NC_TargID(i) && ~NC_RspKey(i)
            NC_miss(s) = NC_miss(s) + 1;
        elseif ~NC_TargID(i) && NC_RspKey(i)
            NC_falseAlarm(s) = NC_falseAlarm(s) + 1;
        elseif ~NC_TargID(i) && ~NC_RspKey(i)
            NC_corReject(s) = NC_corReject(s) + 1;
        end
    end    
 
    NC_H(s) = NC_hit(s)/diff(s);
    if ~NC_H(s)
        NC_H(s) = 0.01;
    end
    NC_F(s) = NC_falseAlarm(s)/same(s);
    if ~NC_F(s)
        NC_F(s) = 0.01;
    end
    NC_dp(s) = norminv(NC_H(s)) - norminv(NC_F(s));
  
    % FOR CROWDING CONDITIONS %%%%%%%
    % KEYPRESS ONLY %%%%
    C_hit_key(s) = 0; C_miss_key(s) = 0; C_falseAlarm_key(s) = 0; C_corReject_key(s) = 0;
    for i = 1:length(C_TargID_key)
        if C_TargID_key(i) && C_RspKey_key(i)
            C_hit_key(s) = C_hit_key(s) + 1;
        elseif C_TargID_key(i) && ~C_RspKey_key(i)
            C_miss_key(s) = C_miss_key(s) + 1;
        elseif ~C_TargID_key(i) && C_RspKey_key(i)
            C_falseAlarm_key(s) = C_falseAlarm_key(s) + 1;
        elseif ~C_TargID_key(i) && ~C_RspKey_key(i)
            C_corReject_key(s) = C_corReject_key(s) + 1;
        end
    end
    
    C_H_key(s) = C_hit_key(s)/diff(s);
    if ~C_H_key(s)
        C_H_key(s) = 0.01;
    end
    C_F_key(s) = C_falseAlarm_key(s)/same(s);
    if ~C_F_key(s)
        C_F_key(s) = 0.01;
    end
    C_dp_key(s) = norminv(C_H_key(s)) - norminv(C_F_key(s));
    
    % GRASP %%%%
    C_hit(s) = 0; C_miss(s) = 0; C_falseAlarm(s) = 0; C_corReject(s) = 0;
    for i = 1:length(C_TargID)
        if C_TargID(i) && C_RspKey(i)
            C_hit(s) = C_hit(s) + 1;
        elseif C_TargID(i) && ~C_RspKey(i)
            C_miss(s) = C_miss(s) + 1;
        elseif ~C_TargID(i) && C_RspKey(i)
            C_falseAlarm(s) = C_falseAlarm(s) + 1;
        elseif ~C_TargID(i) && ~C_RspKey(i)
            C_corReject(s) = C_corReject(s) + 1;
        end
    end   
    
    C_H(s) = C_hit(s)/diff(s);
    if ~C_H(s)
        C_H(s) = 0.01;
    end
    C_F(s) = C_falseAlarm(s)/same(s);
    if ~C_F(s)
        C_F(s) = 0.01;
    end
    C_dp(s) = norminv(C_H(s)) - norminv(C_F(s));
    
end

% COMPARE GRASP/KEY D-PRIMES
% NC_dp_key, NC_dp, C_dp_key, C_dp

NC_diff = NC_dp - NC_dp_key;
C_diff = C_dp - C_dp_key;

%% PLOT THRESHOLDS
close all
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
x1 = thresholds(:,1)';
% x1(length(thresholds) + 1) = mean(thresholds(:,1));
x2 = thresholds(:,2)';
% x2(length(thresholds) + 1) = mean(thresholds(:,2));

h = bar([x1;x2]);
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
% h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northwest');
set(gca,gcaOpts{:})
title('Orientation thresholds')
xlabel('Crowding condition')
ylabel('Threshold (deg)')
ylim([0 6])

%% PLOT HITS etc KEY
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:8,'XTickLabel',{'NC hit' 'NC miss' 'NC FA' 'NC CR' 'C hit' 'C miss' 'C FA' 'C CR'},'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([NC_hit_key; NC_miss_key; NC_falseAlarm_key; NC_corReject_key; ...
    C_hit_key; C_miss_key; C_falseAlarm_key; C_corReject_key]); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
% h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('Orientation sensitivity (d'')')
xlabel('Condition')
ylabel('d-prime')
ylim([0 48])

%% PLOT HITS etc GRASP
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:8,'XTickLabel',{'NC hit' 'NC miss' 'NC FA' 'NC CR' 'C hit' 'C miss' 'C FA' 'C CR'},'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([NC_hit; NC_miss; NC_falseAlarm; NC_corReject; ...
    C_hit; C_miss; C_falseAlarm; C_corReject]); 
% h(1).FaceColor = [135,205,215]/255;
% h(2).FaceColor = [240,59,37]/255;
% h(3).FaceColor = [250,185,219]/255;
% h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('GRASP')
xlabel('Condition')
ylabel('d-prime')
ylim([0 50])


%% PLOT D PRIMES / WITH GRASP
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([NC_dp_key; NC_dp; C_dp_key; C_dp]); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
% h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('Orientation sensitivity (d'')')
xlabel('Condition')
ylabel('d-prime')
ylim([-1 3.5])

%% PLOT D PRIME COMPARISON, GRASP VS NO GRASP
cndNames = {'NC','C'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([NC_diff; C_diff]); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
% h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('d'' diff (grasping - key only)')
xlabel('Condition')
ylabel('d'' difference')
ylim([0 2])

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
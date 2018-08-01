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
% 4 Ryan - 1.25x multiplier

% calculate reaction time after s2 onset + initiation latency for movement

%%

xy1 = data.xy1;
xy2 = data.xy2;

xy1_data1 = data.xy1_data1;
xy1_data2 = data.xy1_data2;
xy2_data1 = data.xy2_data1;
xy2_data2 = data.xy2_data2;

% figure
% hold on
% plot(xy1_data1,'color','r');
% plot(xy1_data2,'color','m');
% plot(xy2_data1,'color','b');
% plot(xy2_data2,'color','c');

figure
hold on
plot(xy1_data1,xy1_data2*-1,'o','color','r')
plot(xy2_data1,xy2_data2*-1,'o','color','b')

%%

radius = sqrt((xy2(1)-xy1(1))^2 + (xy2(2)-xy1(2))^2); % distance between two fingers
x_dist = xy2(1) - xy1(1); % x dist between two trackers 
angle = abs(rad2deg(asin(x_dist/radius))); % angle between two fingers
if baseOri(trials) % IF HORIZONTAL
    angle = 90 - angle;
end

%%
clear all
close all

% WHICH DATASET?
dataset = 1; % 0 = main exp no grasp, 1 = main exp with grasp

% INCLUDE CATCH?
incCatch = 0;

% parentDir = '~/Bethany/paclab';
parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));
dataDir = sprintf('%s/Subject_folders',parentDir);
if ~dataset
    dataDir = sprintf('%s/Pilot_NoGrasp',dataDir);
end
subj = sort(strsplit(ls(dataDir)));
% subj = subj(7:end);
subj = subj(end);
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
diff = zeros(1,nSubj); same = zeros(1,nSubj); 
nInvalidKey = zeros(nSubj,2); percAccGrasp = zeros(nSubj,2);

RH_NC_TargID_key = [];
RH_NC_RspKey_key = [];
RH_C_TargID_key = [];
RH_C_RspKey_key = [];

if ~incCatch
    catch_NC_hit_key = zeros(1,nSubj); catch_NC_miss_key = zeros(1,nSubj); catch_NC_falseAlarm_key = zeros(1,nSubj); 
    catch_NC_corReject_key = zeros(1,nSubj);
    catch_NC_H_key = zeros(1,nSubj); catch_NC_F_key = zeros(1,nSubj); catch_NC_dp_key = zeros(1,nSubj);
    catch_C_hit_key = zeros(1,nSubj); catch_C_miss_key = zeros(1,nSubj); catch_C_falseAlarm_key = zeros(1,nSubj); 
    catch_C_corReject_key = zeros(1,nSubj);
    catch_C_H_key = zeros(1,nSubj); catch_C_F_key = zeros(1,nSubj); catch_C_dp_key = zeros(1,nSubj);
    catch_NC_hit = zeros(1,nSubj); catch_NC_miss = zeros(1,nSubj); catch_NC_falseAlarm = zeros(1,nSubj); 
    catch_NC_corReject = zeros(1,nSubj);
    catch_NC_H = zeros(1,nSubj); catch_NC_F = zeros(1,nSubj); catch_NC_dp = zeros(1,nSubj);
    catch_C_hit = zeros(1,nSubj); catch_C_miss = zeros(1,nSubj); catch_C_falseAlarm = zeros(1,nSubj); 
    catch_C_corReject = zeros(1,nSubj);
    catch_C_H = zeros(1,nSubj); catch_C_F = zeros(1,nSubj); catch_C_dp = zeros(1,nSubj);
    catch_diff = zeros(1,nSubj); catch_same = zeros(1,nSubj); 
end

for s = 1:nSubj
    NC_TargID = []; NC_RspKey = []; C_TargID = []; C_RspKey = [];
    NC_TargID_key = []; NC_RspKey_key = []; C_TargID_key = []; C_RspKey_key = [];
    if ~incCatch
        catch_NC_TargID = []; catch_NC_RspKey = []; catch_C_TargID = []; catch_C_RspKey = [];
        catch_NC_TargID_key = []; catch_NC_RspKey_key = []; catch_C_TargID_key = []; catch_C_RspKey_key = [];
    end
    
    curSubj = subj{s};
    subjNum = str2num(curSubj(2:3));
    subjDir = sprintf('%s/%s',dataDir,curSubj);
    blockFiles = mySubFiles(subjDir,'all',13);
%     blockFiles = mySubFiles(subjDir,'.mat',12);
    txtFiles = mySubFiles(subjDir,'.txt',13);
    
    dataThresh = load(blockFiles{1});
    thresholds(s,:) = dataThresh.thresholds;
    
    idx = 1; grIdx = 1;blockIdx = 1;
    for b = 2:length(blockFiles)
        data = load(blockFiles{b});
        targID = data.targID;
        baseOri = data.baseOri;
        rspKey = data.rspKey;
        allCndList = data.cndList;
        catchIdx = data.catchIdx;
        nonCatchIdx = data.nonCatchIdx;
        nActual = data.trialNumber;
        nCatch = data.nCatch;
        
        horizIdx = find(baseOri);
        vertIdx = find(~baseOri);
        
        acc = data.acc;
        
        % RT ANALYSIS (for all trials)
        S1_onset_time = data.S1_onset_time;
        S2_onset_time = data.S2_onset_time;
        rt = data.rt;

        if ~data.grasping
           nInvalidKey(s,idx) = length(find(rt<(.1 + .3 + .1 + .1))); 
           idx = idx + 1;
        else
           SOT_data = data.SOT_data;
           percAccGrasp(s,grIdx) = sum(data.acc_grasp)/length(data.acc_grasp);
%            timeElapsed = data.timeElapsed - (.1 + .3 + .1); % ONLY IF TIME ELAPSED EXISTS
%            initLatency = SOT_data(end) - S1_onset_time 
           grIdx = grIdx + 1;
        end
        
        % TRAJECTORY ANALYSIS
        if data.grasping 
            clear trialIdx
            curBlock = dlmread(txtFiles{blockIdx});
            xy1_data1{blockIdx,s} = curBlock(:,2);
            xy1_data2{blockIdx,s} = curBlock(:,3);
            xy2_data1{blockIdx,s} = curBlock(:,4);
            xy2_data2{blockIdx,s} = curBlock(:,5);
            
            idx = 1;
            for x = 1:length(xy1_data1{blockIdx,s})-1
                if abs(xy1_data1{blockIdx,s}(x+1) - xy1_data1{blockIdx,s}(x)) > 200
                    trialIdx(idx) = x; % gets the last sample of the trial
                    idx = idx + 1;
                end
            end
            trialIdx(idx) = length(xy1_data1{blockIdx,s});
            
            final_xy1_d1(:,blockIdx,s) = xy1_data1{blockIdx,s}(trialIdx);
            final_xy1_d2(:,blockIdx,s) = xy1_data2{blockIdx,s}(trialIdx);
            final_xy2_d1(:,blockIdx,s) = xy2_data1{blockIdx,s}(trialIdx);
            final_xy2_d2(:,blockIdx,s) = xy2_data2{blockIdx,s}(trialIdx);
            
%             actual_idx = find(final_xy1_d1(:,blockIdx,s) > 700);
            actual_xy1_d1(:,blockIdx,s) = final_xy1_d1(nonCatchIdx,blockIdx,s); % len = 48 bc 48 real trials per 64-trial block
            actual_xy1_d2(:,blockIdx,s) = final_xy1_d2(nonCatchIdx,blockIdx,s);
            actual_xy2_d1(:,blockIdx,s) = final_xy2_d1(nonCatchIdx,blockIdx,s);
            actual_xy2_d2(:,blockIdx,s) = final_xy2_d2(nonCatchIdx,blockIdx,s);
            
            trialEndpoint(:,blockIdx,s) = trialIdx;
            blockIdx = blockIdx + 1;
        end
        
        allTrials = 1:length(targID); % 128 trials, per block
        
        % targID & allCndList = complete 64 trials
        % nonCatchIdx = index of all non-catch in the full 64
        
        all_NC_idx = sort(find(allCndList == 5 | allCndList == 6 | allCndList == 7 | allCndList == 8)); % length = 48 / non crowd
        all_C_idx = sort(find(allCndList == 1 | allCndList == 2 | allCndList == 3 | allCndList == 4 ...
             | allCndList == 9 | allCndList == 10 | allCndList == 11 | allCndList == 12));
        
        RH_NC_idx = intersect(all_NC_idx,nonCatchIdx);
        RH_C_idx = intersect(all_C_idx,nonCatchIdx);
        LH_NC_idx = intersect(all_NC_idx,catchIdx);
        LH_C_idx = intersect(all_C_idx,catchIdx);
        
        x = targID(RH_NC_idx); % NONCROWDING TARGET ID; length 24 (post-ryan)
        y = rspKey(RH_NC_idx); % NONCROWDING RESP KEY; length 24 (post-ryan)
        z = targID(RH_C_idx); % CROWDING TARGET ID; length 24 (post-ryan)
        a = rspKey(RH_C_idx); % CROWDING RESP KEY; length 24 (post-ryan)
        
        r = targID(LH_NC_idx); % NONCROWDING TARGET ID; length 24 (post-ryan)
        t = rspKey(LH_NC_idx); % NONCROWDING RESP KEY; length 24 (post-ryan)
        u = targID(LH_C_idx); % CROWDING TARGET ID; length 24 (post-ryan)
        v = rspKey(LH_C_idx); % CROWDING RESP KEY; length 24 (post-ryan)
        
        if ~data.grasping
            RH_NC_TargID_key = [RH_NC_TargID_key targID(RH_NC_idx)];
            RH_NC_RspKey_key = [RH_NC_RspKey_key rspKey(RH_NC_idx)];
            RH_C_TargID_key = [RH_C_TargID_key targID(RH_C_idx)];
            RH_C_RspKey_key = [RH_C_RspKey_key rspKey(RH_C_idx)];
            
            NC_TargID_key = cat(2,NC_TargID_key,x); % len = 96, because 48 per block * 2 blocks
            NC_RspKey_key = cat(2,NC_RspKey_key,y); % len = 96, because 48 per block * 2 blocks
            C_TargID_key = cat(2,C_TargID_key,z); % len = 96, because 48 per block * 2 blocks
            C_RspKey_key = cat(2,C_RspKey_key,a); % len = 96, because 48 per block * 2 blocks  
        else
            NC_TargID = cat(2,NC_TargID,x); % len = 96, because 48 per block * 2 blocks
            NC_RspKey = cat(2,NC_RspKey,y); % len = 96, because 48 per block * 2 blocks
            C_TargID = cat(2,C_TargID,z); % len = 96, because 48 per block * 2 blocks
            C_RspKey = cat(2,C_RspKey,a); % len = 96, because 48 per block * 2 blocks
        end
    end % END OF BLOCKWISE ANALYSIS
    
%     diff(s) = sum(NC_TargID); % = 48, because 48 total DIFF trials per task type (grasp/key) per subj; 24 * 2
%     same(s) = length(NC_TargID) - sum(NC_TargID); % = 48, same reasoning ^
%     if ~incCatch
%         catch_diff(s) = sum(catch_NC_TargID);% = 16, because 16 total DIFF trials per task type per subj; 8 * 2
%         catch_same(s) = length(catch_NC_TargID) - sum(catch_NC_TargID); % = 16, same reasoning ^
%     end
%     
%     % FOR NO CROWDING CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % KEYPRESS ONLY %%%%
%     NC_hit_key(s) = 0; NC_miss_key(s) = 0; NC_falseAlarm_key(s) = 0; NC_corReject_key(s) = 0;
%     for i = 1:length(NC_TargID_key)
%         if NC_TargID_key(i) && NC_RspKey_key(i)
%             NC_hit_key(s) = NC_hit_key(s) + 1;
%         elseif NC_TargID_key(i) && ~NC_RspKey_key(i)
%             NC_miss_key(s) = NC_miss_key(s) + 1;
%         elseif ~NC_TargID_key(i) && NC_RspKey_key(i)
%             NC_falseAlarm_key(s) = NC_falseAlarm_key(s) + 1;
%         elseif ~NC_TargID_key(i) && ~NC_RspKey_key(i)
%             NC_corReject_key(s) = NC_corReject_key(s) + 1;
%         end
%     end
%     
%     NC_H_key(s) = NC_hit_key(s)/diff(s);
%     if ~NC_H_key(s)
%         NC_H_key(s) = 0.01;
%     end
%     NC_F_key(s) = NC_falseAlarm_key(s)/same(s);
%     if ~NC_F_key(s)
%         NC_F_key(s) = 0.01;
%     end
%     NC_dp_key(s) = norminv(NC_H_key(s)) - norminv(NC_F_key(s));
%     
%     % GRASP %%%%
%     NC_hit(s) = 0; NC_miss(s) = 0; NC_falseAlarm(s) = 0; NC_corReject(s) = 0;
%     for i = 1:length(NC_TargID)
%         if NC_TargID(i) && NC_RspKey(i)
%             NC_hit(s) = NC_hit(s) + 1;
%         elseif NC_TargID(i) && ~NC_RspKey(i)
%             NC_miss(s) = NC_miss(s) + 1;
%         elseif ~NC_TargID(i) && NC_RspKey(i)
%             NC_falseAlarm(s) = NC_falseAlarm(s) + 1;
%         elseif ~NC_TargID(i) && ~NC_RspKey(i)
%             NC_corReject(s) = NC_corReject(s) + 1;
%         end
%     end    
%  
%     NC_H(s) = NC_hit(s)/diff(s);
%     if ~NC_H(s)
%         NC_H(s) = 0.01;
%     end
%     NC_F(s) = NC_falseAlarm(s)/same(s);
%     if ~NC_F(s)
%         NC_F(s) = 0.01;
%     end
%     NC_dp(s) = norminv(NC_H(s)) - norminv(NC_F(s));
%   
%     % FOR CROWDING CONDITIONS %%%%%%%
%     % KEYPRESS ONLY %%%%
%     C_hit_key(s) = 0; C_miss_key(s) = 0; C_falseAlarm_key(s) = 0; C_corReject_key(s) = 0;
%     for i = 1:length(C_TargID_key)
%         if C_TargID_key(i) && C_RspKey_key(i)
%             C_hit_key(s) = C_hit_key(s) + 1;
%         elseif C_TargID_key(i) && ~C_RspKey_key(i)
%             C_miss_key(s) = C_miss_key(s) + 1;
%         elseif ~C_TargID_key(i) && C_RspKey_key(i)
%             C_falseAlarm_key(s) = C_falseAlarm_key(s) + 1;
%         elseif ~C_TargID_key(i) && ~C_RspKey_key(i)
%             C_corReject_key(s) = C_corReject_key(s) + 1;
%         end
%     end
%     
%     C_H_key(s) = C_hit_key(s)/diff(s);
%     if ~C_H_key(s)
%         C_H_key(s) = 0.01;
%     end
%     C_F_key(s) = C_falseAlarm_key(s)/same(s);
%     if ~C_F_key(s)
%         C_F_key(s) = 0.01;
%     end
%     C_dp_key(s) = norminv(C_H_key(s)) - norminv(C_F_key(s));
%     
%     % GRASP %%%%
%     C_hit(s) = 0; C_miss(s) = 0; C_falseAlarm(s) = 0; C_corReject(s) = 0;
%     for i = 1:length(C_TargID)
%         if C_TargID(i) && C_RspKey(i)
%             C_hit(s) = C_hit(s) + 1;
%         elseif C_TargID(i) && ~C_RspKey(i)
%             C_miss(s) = C_miss(s) + 1;
%         elseif ~C_TargID(i) && C_RspKey(i)
%             C_falseAlarm(s) = C_falseAlarm(s) + 1;
%         elseif ~C_TargID(i) && ~C_RspKey(i)
%             C_corReject(s) = C_corReject(s) + 1;
%         end
%     end   
%     
%     C_H(s) = C_hit(s)/diff(s);
%     if ~C_H(s)
%         C_H(s) = 0.01;
%     end
%     C_F(s) = C_falseAlarm(s)/same(s);
%     if ~C_F(s)
%         C_F(s) = 0.01;
%     end
%     C_dp(s) = norminv(C_H(s)) - norminv(C_F(s));
%     
%     
%     % CATCH TRIALS! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % FOR NO CROWDING CONDITIONS %%%%%%%%%
%     % KEYPRESS ONLY %%%%
%     if ~incCatch
%         catch_NC_hit_key(s) = 0; catch_NC_miss_key(s) = 0; catch_NC_falseAlarm_key(s) = 0; catch_NC_corReject_key(s) = 0;
%         for i = 1:length(catch_NC_TargID_key)
%             if catch_NC_TargID_key(i) && catch_NC_RspKey_key(i)
%                 catch_NC_hit_key(s) = catch_NC_hit_key(s) + 1;
%             elseif catch_NC_TargID_key(i) && ~catch_NC_RspKey_key(i)
%                 catch_NC_miss_key(s) = catch_NC_miss_key(s) + 1;
%             elseif ~catch_NC_TargID_key(i) && catch_NC_RspKey_key(i)
%                 catch_NC_falseAlarm_key(s) = catch_NC_falseAlarm_key(s) + 1;
%             elseif ~catch_NC_TargID_key(i) && ~catch_NC_RspKey_key(i)
%                 catch_NC_corReject_key(s) = catch_NC_corReject_key(s) + 1;
%             end
%         end
% 
%         catch_NC_H_key(s) = catch_NC_hit_key(s)/catch_diff(s);
%         if ~catch_NC_H_key(s)
%             catch_NC_H_key(s) = 0.01;
%         end
%         catch_NC_F_key(s) = catch_NC_falseAlarm_key(s)/catch_same(s);
%         if ~catch_NC_F_key(s)
%             catch_NC_F_key(s) = 0.01;
%         end
%         catch_NC_dp_key(s) = norminv(catch_NC_H_key(s)) - norminv(catch_NC_F_key(s));
% 
%         % GRASP %%%%
%         catch_NC_hit(s) = 0; catch_NC_miss(s) = 0; catch_NC_falseAlarm(s) = 0; catch_NC_corReject(s) = 0;
%         for i = 1:length(catch_NC_TargID)
%             if catch_NC_TargID(i) && catch_NC_RspKey(i)
%                 catch_NC_hit(s) = catch_NC_hit(s) + 1;
%             elseif catch_NC_TargID(i) && ~catch_NC_RspKey(i)
%                 catch_NC_miss(s) = catch_NC_miss(s) + 1;
%             elseif ~catch_NC_TargID(i) && catch_NC_RspKey(i)
%                 catch_NC_falseAlarm(s) = catch_NC_falseAlarm(s) + 1;
%             elseif ~catch_NC_TargID(i) && ~catch_NC_RspKey(i)
%                 catch_NC_corReject(s) = catch_NC_corReject(s) + 1;
%             end
%         end    
% 
%         catch_NC_H(s) = catch_NC_hit(s)/catch_diff(s);
%         if ~catch_NC_H(s)
%             catch_NC_H(s) = 0.01;
%         end
%         catch_NC_F(s) = catch_NC_falseAlarm(s)/catch_same(s);
%         if ~catch_NC_F(s)
%             catch_NC_F(s) = 0.01;
%         end
%         catch_NC_dp(s) = norminv(catch_NC_H(s)) - norminv(catch_NC_F(s));
% 
%         % FOR CROWDING CONDITIONS %%%%%%%
%         % KEYPRESS ONLY %%%%
%         catch_C_hit_key(s) = 0; catch_C_miss_key(s) = 0; catch_C_falseAlarm_key(s) = 0; catch_C_corReject_key(s) = 0;
%         for i = 1:length(catch_C_TargID_key)
%             if catch_C_TargID_key(i) && catch_C_RspKey_key(i)
%                 catch_C_hit_key(s) = catch_C_hit_key(s) + 1;
%             elseif catch_C_TargID_key(i) && ~catch_C_RspKey_key(i)
%                 catch_C_miss_key(s) = catch_C_miss_key(s) + 1;
%             elseif ~catch_C_TargID_key(i) && catch_C_RspKey_key(i)
%                 catch_C_falseAlarm_key(s) = catch_C_falseAlarm_key(s) + 1;
%             elseif ~catch_C_TargID_key(i) && ~catch_C_RspKey_key(i)
%                 catch_C_corReject_key(s) = catch_C_corReject_key(s) + 1;
%             end
%         end
% 
%         catch_C_H_key(s) = catch_C_hit_key(s)/catch_diff(s);
%         if ~catch_C_H_key(s)
%             catch_C_H_key(s) = 0.01;
%         end
%         catch_C_F_key(s) = catch_C_falseAlarm_key(s)/catch_same(s);
%         if ~catch_C_F_key(s)
%             catch_C_F_key(s) = 0.01;
%         end
%         catch_C_dp_key(s) = norminv(catch_C_H_key(s)) - norminv(catch_C_F_key(s));
% 
%         % GRASP %%%%
%         catch_C_hit(s) = 0; catch_C_miss(s) = 0; catch_C_falseAlarm(s) = 0; catch_C_corReject(s) = 0;
%         for i = 1:length(catch_C_TargID)
%             if catch_C_TargID(i) && catch_C_RspKey(i)
%                 catch_C_hit(s) = catch_C_hit(s) + 1;
%             elseif catch_C_TargID(i) && ~catch_C_RspKey(i)
%                 catch_C_miss(s) = catch_C_miss(s) + 1;
%             elseif ~catch_C_TargID(i) && catch_C_RspKey(i)
%                 catch_C_falseAlarm(s) = catch_C_falseAlarm(s) + 1;
%             elseif ~catch_C_TargID(i) && ~catch_C_RspKey(i)
%                 catch_C_corReject(s) = catch_C_corReject(s) + 1;
%             end
%         end   
% 
%         catch_C_H(s) = catch_C_hit(s)/catch_diff(s);
%         if ~catch_C_H(s)
%             catch_C_H(s) = 0.01;
%         end
%         catch_C_F(s) = catch_C_falseAlarm(s)/catch_same(s);
%         if ~catch_C_F(s)
%             catch_C_F(s) = 0.01;
%         end
%         catch_C_dp(s) = norminv(catch_C_H(s)) - norminv(catch_C_F(s));
%     end
end

% COMPARE GRASP/KEY D-PRIMES
% NC_dp_key, NC_dp, C_dp_key, C_dp

% NC_diff = NC_dp - NC_dp_key;
% C_diff = C_dp - C_dp_key;
% if ~incCatch
%     catch_NC_diff = catch_NC_dp - catch_NC_dp_key;
%     catch_C_diff = catch_C_dp - catch_C_dp_key;
% end

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
h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northwest');
set(gca,gcaOpts{:})
title('Orientation thresholds')
xlabel('Crowding condition')
ylabel('Threshold (deg)')
ylim([0 6])

%% PLOT TRAJECTORIES
gcaOpts = {'box','off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
for s = 1:nSubj
    subplot(1,nSubj,s)
    hold on
    for b = 1:4
        plot(xy1_data1{b,s}(1:trialEndpoint(1,b,s)),xy1_data2{b,s}(1:trialEndpoint(1,b,s))*-1,'b')
        plot(xy2_data1{b,s}(1:trialEndpoint(1,b,s)),xy2_data2{b,s}(1:trialEndpoint(1,b,s))*-1,'r')
        for u = 1:length(trialEndpoint(:,b,s))-1
            plot(xy1_data1{b,s}(trialEndpoint(u,b,s)+1:trialEndpoint(u+1,b,s)),xy1_data2{b,s}(trialEndpoint(u,b,s)+...
                1:trialEndpoint(u+1,b,s))*-1,'color','b')
            plot(xy2_data1{b,s}(trialEndpoint(u,b,s)+1:trialEndpoint(u+1,b,s)),xy2_data2{b,s}(trialEndpoint(u,b,s)+...
                1:trialEndpoint(u+1,b,s))*-1,'color','r')
        end
    end
end
set(gca,gcaOpts{:})
ylim([-1250 -300]);
xlim([-200 1400]);
title('Trajectories for all trials')

% SPECIFIC COMBINATIONS
% crowded + horizontal
% crowded + vertical
% noncrowded + horizontal
% noncrowded + vertical

subsetIdx = [];
subsetEndpoint = trialEndpoint(subsetIdx,:,:);

C_trialEndpoint = trialEndpoint(C_idx,:,:);
figure
for s = 1:nSubj
    subplot(1,nSubj,s)
    hold on
    for b = 1:4
        plot(xy1_data1{b,s}(1:C_trialEndpoint(1,b,s)),xy1_data2{b,s}(1:C_trialEndpoint(1,b,s))*-1,'b')
        plot(xy2_data1{b,s}(1:C_trialEndpoint(1,b,s)),xy2_data2{b,s}(1:C_trialEndpoint(1,b,s))*-1,'r')
        for u = 1:length(C_trialEndpoint(:,b,s))-1
            plot(xy1_data1{b,s}(C_trialEndpoint(u,b,s)+1:C_trialEndpoint(u+1,b,s)),xy1_data2{b,s}(C_trialEndpoint(u,b,s)+...
                1:C_trialEndpoint(u+1,b,s))*-1,'color','b')
            plot(xy2_data1{b,s}(C_trialEndpoint(u,b,s)+1:C_trialEndpoint(u+1,b,s)),xy2_data2{b,s}(C_trialEndpoint(u,b,s)+...
                1:C_trialEndpoint(u+1,b,s))*-1,'color','r')
        end
    end
end
set(gca,gcaOpts{:})
ylim([-1250 -300]);
xlim([-200 1400]);
title('Trajectories for crowding trials')

NC_trialEndpoint = trialEndpoint(NC_idx,:,:);
figure
for s = 1:nSubj
    subplot(1,nSubj,s)
    hold on
    for b = 1:4
        plot(xy1_data1{b,s}(1:NC_trialEndpoint(1,b,s)),xy1_data2{b,s}(1:NC_trialEndpoint(1,b,s))*-1,'b')
        plot(xy2_data1{b,s}(1:NC_trialEndpoint(1,b,s)),xy2_data2{b,s}(1:NC_trialEndpoint(1,b,s))*-1,'r')
        for u = 1:length(NC_trialEndpoint(:,b,s))-1
            plot(xy1_data1{b,s}(NC_trialEndpoint(u,b,s)+1:NC_trialEndpoint(u+1,b,s)),xy1_data2{b,s}(NC_trialEndpoint(u,b,s)+...
                1:NC_trialEndpoint(u+1,b,s))*-1,'color','b')
            plot(xy2_data1{b,s}(NC_trialEndpoint(u,b,s)+1:NC_trialEndpoint(u+1,b,s)),xy2_data2{b,s}(NC_trialEndpoint(u,b,s)+...
                1:NC_trialEndpoint(u+1,b,s))*-1,'color','r')
        end
    end
end
set(gca,gcaOpts{:})
ylim([-1250 -300]);
xlim([-200 1400]);
title('Trajectories for non-crowding trials')


%% PLOT HITS etc KEY
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:8,'XTickLabel',{'NC hit' 'NC miss' 'NC FA' 'NC CR' 'C hit' 'C miss' 'C FA' 'C CR'},'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([catch_NC_hit_key; catch_NC_miss_key; catch_NC_falseAlarm_key; catch_NC_corReject_key; ...
    catch_C_hit_key; catch_C_miss_key; catch_C_falseAlarm_key; catch_C_corReject_key]); 
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
h = bar([catch_NC_hit; catch_NC_miss; catch_NC_falseAlarm; catch_NC_corReject; ...
    catch_C_hit; catch_C_miss; catch_C_falseAlarm; catch_C_corReject]); 
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

%% PLOT D PRIMES / WITH GRASP / FOR NON-CATCH
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([NC_dp_key; NC_dp; C_dp_key; C_dp]); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('Orientation sensitivity (d'')')
xlabel('Condition')
ylabel('d-prime')
ylim([-1 3.5])


%% PLOT D PRIMES / WITH GRASP / FOR CATCH
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([catch_NC_dp_key; catch_NC_dp; catch_C_dp_key; catch_C_dp]); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('Orientation sensitivity (d'') for CATCH TRIALS')
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
h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('d'' diff (grasping - key only)')
xlabel('Condition')
ylabel('d'' difference')
ylim([-2 2])

%% PLOT D PRIME COMPARISON, GRASP VS NO GRASP / CATCH TRIALS
cndNames = {'NC','C'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([catch_NC_diff; catch_C_diff]); 
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;
legend({subj{:}},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('d'' diff (grasping - key only) for CATCH TRIALS')
xlabel('Condition')
ylabel('d'' difference')
ylim([-2 2])

%% PLOT D PRIMES / NO GRASP
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
h = bar([catch_NC_dp_key, mean(catch_NC_dp_key); catch_C_dp_key, mean(catch_C_dp_key)]); 
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
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
clear all
close all

% WHICH DATASET?
dataset = 1; % 0 = main exp no grasp, 1 = main exp with grasp

% INCLUDE CATCH?
incCatch = 0;

parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));
dataDir = sprintf('%s/Subject_folders',parentDir);
if ~dataset
    dataDir = sprintf('%s/Pilot_NoGrasp',dataDir);
end
subj = sort(strsplit(ls(dataDir)));
subj = subj(7:end);
% subj = subj(end);
nSubj = length(subj);

% PREALLOCATION
nInvalidKey = zeros(nSubj,2); percAccGrasp = zeros(nSubj,2);
grasping = nan(8,nSubj);
horizIdx = nan(64,8,nSubj); vertIdx = nan(64,8,nSubj);
% cndList = nan(128,8,nSubj);
LH_idx = cell(8,nSubj);RH_idx = cell(8,nSubj);
H_idx = cell(8,nSubj);V_idx = cell(8,nSubj);
NC_idx = cell(8,nSubj);C_idx = cell(8,nSubj);
LH_H_NC_idx = cell(8,nSubj);LH_H_C_idx = cell(8,nSubj);LH_V_NC_idx = cell(8,nSubj);LH_V_C_idx = cell(8,nSubj);
RH_H_NC_idx = cell(8,nSubj);RH_H_C_idx = cell(8,nSubj);RH_V_NC_idx = cell(8,nSubj);RH_V_C_idx = cell(8,nSubj);

for s = 1:nSubj
    RH_NC_Ans_key = [];RH_NC_Rsp_key = [];RH_C_Ans_key = [];RH_C_Rsp_key = [];
    LH_NC_Ans_key = [];LH_NC_Rsp_key = [];LH_C_Ans_key = [];LH_C_Rsp_key = [];
    RH_NC_Ans_grp = [];RH_NC_Rsp_grp = [];RH_C_Ans_grp = [];RH_C_Rsp_grp = [];
    LH_NC_Ans_grp = [];LH_NC_Rsp_grp = [];LH_C_Ans_grp = [];LH_C_Rsp_grp = [];
    
    RH_NC_Ans_key_old = [];
    
    curSubj = subj{s};
    subjNum = str2num(curSubj(2:3));
    subjDir = sprintf('%s/%s',dataDir,curSubj);
    blockFiles = mySubFiles(subjDir,'all',13);
%     blockFiles = mySubFiles(subjDir,'.mat',12);
    if s >= 4
        txtFiles = mySubFiles(subjDir,'.txt',13);
    end
    
    dataThresh = load(blockFiles{1});
    thresholds(s,:) = dataThresh.thresholds;
    
    idx = 1; grIdx = 1;blockIdx = 1;
    for b = 2:length(blockFiles)
        data = load(blockFiles{b});
        grasping(b-1,s) = data.grasping;
        targID = data.targID;
        baseOri = data.baseOri;
        flanker = data.flanker;
        rspKey = data.rspKey;
        cndList = data.cndList; % (1:length(data.cndList),b-1,s) 
        
        nRH = data.trialNumber;
        nLH = data.nCatch;
        acc = data.acc;
        
        
        allTrials = 1:length(targID); % 128 trials, per block

        % INDICES
        LH_idx{b-1,s} = sort(data.catchIdx);
        RH_idx{b-1,s} = sort(data.nonCatchIdx);
        
        H_idx{b-1,s} = sort(find(baseOri));
        V_idx{b-1,s} = sort(find(~baseOri));
        
        NC_idx{b-1,s} = sort(find(~flanker));
        C_idx{b-1,s} = sort(find(flanker));
        
        LH_H_NC_idx{b-1,s} = mintersect(LH_idx{b-1,s}, H_idx{b-1,s}, NC_idx{b-1,s});
        LH_H_C_idx{b-1,s} = mintersect(LH_idx{b-1,s}, H_idx{b-1,s}, C_idx{b-1,s});
        LH_V_NC_idx{b-1,s} = mintersect(LH_idx{b-1,s}, V_idx{b-1,s}, NC_idx{b-1,s});
        LH_V_C_idx{b-1,s} = mintersect(LH_idx{b-1,s}, V_idx{b-1,s}, C_idx{b-1,s});
        RH_H_NC_idx{b-1,s} = mintersect(RH_idx{b-1,s}, H_idx{b-1,s}, NC_idx{b-1,s});
        RH_H_C_idx{b-1,s} = mintersect(RH_idx{b-1,s}, H_idx{b-1,s}, C_idx{b-1,s});
        RH_V_NC_idx{b-1,s} = mintersect(RH_idx{b-1,s}, V_idx{b-1,s}, NC_idx{b-1,s});
        RH_V_C_idx{b-1,s} = mintersect(RH_idx{b-1,s}, V_idx{b-1,s}, C_idx{b-1,s});
        
        H_NC_idx{b-1,s} = intersect(H_idx{b-1,s}, NC_idx{b-1,s});
        H_C_idx{b-1,s} = intersect(H_idx{b-1,s}, C_idx{b-1,s});
        V_NC_idx{b-1,s} = intersect(V_idx{b-1,s}, NC_idx{b-1,s});
        V_C_idx{b-1,s} = intersect(V_idx{b-1,s}, C_idx{b-1,s});
        
        if ~data.grasping
            RH_NC_Ans_key = [RH_NC_Ans_key targID([RH_H_NC_idx{b-1,s},RH_V_NC_idx{b-1,s}])]; % len = 96; 24/block * 4 blocks
            RH_NC_Rsp_key = [RH_NC_Rsp_key rspKey([RH_H_NC_idx{b-1,s},RH_V_NC_idx{b-1,s}])];
            RH_C_Ans_key = [RH_C_Ans_key targID([RH_H_C_idx{b-1,s},RH_V_C_idx{b-1,s}])];
            RH_C_Rsp_key = [RH_C_Rsp_key rspKey([RH_H_C_idx{b-1,s},RH_V_C_idx{b-1,s}])];
            LH_NC_Ans_key = [LH_NC_Ans_key targID([LH_H_NC_idx{b-1,s},LH_V_NC_idx{b-1,s}])];
            LH_NC_Rsp_key = [LH_NC_Rsp_key rspKey([LH_H_NC_idx{b-1,s},LH_V_NC_idx{b-1,s}])];
            LH_C_Ans_key = [LH_C_Ans_key targID([LH_H_C_idx{b-1,s},LH_V_C_idx{b-1,s}])];
            LH_C_Rsp_key = [LH_C_Rsp_key rspKey([LH_H_C_idx{b-1,s},LH_V_C_idx{b-1,s}])];
        else
            RH_NC_Ans_grp = [RH_NC_Ans_grp targID([RH_H_NC_idx{b-1,s},RH_V_NC_idx{b-1,s}])]; % len = 96; 24/block * 4 blocks
            RH_NC_Rsp_grp = [RH_NC_Rsp_grp rspKey([RH_H_NC_idx{b-1,s},RH_V_NC_idx{b-1,s}])];
            RH_C_Ans_grp = [RH_C_Ans_grp targID([RH_H_C_idx{b-1,s},RH_V_C_idx{b-1,s}])];
            RH_C_Rsp_grp = [RH_C_Rsp_grp rspKey([RH_H_C_idx{b-1,s},RH_V_C_idx{b-1,s}])];
            LH_NC_Ans_grp = [LH_NC_Ans_grp targID([LH_H_NC_idx{b-1,s},LH_V_NC_idx{b-1,s}])]; 
            LH_NC_Rsp_grp = [LH_NC_Rsp_grp rspKey([LH_H_NC_idx{b-1,s},LH_V_NC_idx{b-1,s}])];
            LH_C_Ans_grp = [LH_C_Ans_grp targID([LH_H_C_idx{b-1,s},LH_V_C_idx{b-1,s}])];
            LH_C_Rsp_grp = [LH_C_Rsp_grp rspKey([LH_H_C_idx{b-1,s},LH_V_C_idx{b-1,s}])];
        end
        
%         TRAJECTORY ANALYSIS
        if s >= 4 && data.grasping 
            angleInit(:,b-1) = data.angleGraspInitial;
            angleFinal(:,b-1) = data.angleGraspFinal;
            angleDiff(:,b-1) = angleFinal(:,b-1) - angleInit(:,b-1);
%         if data.grasping 
            clear trialIdx
            curBlock = dlmread(txtFiles{blockIdx});
            xy1_data1{blockIdx,s-3} = curBlock(:,2);
            xy1_data2{blockIdx,s-3} = curBlock(:,3);
            xy2_data1{blockIdx,s-3} = curBlock(:,4);
            xy2_data2{blockIdx,s-3} = curBlock(:,5);
%             xy1_data1{blockIdx,s} = curBlock(:,2);
%             xy1_data2{blockIdx,s} = curBlock(:,3);
%             xy2_data1{blockIdx,s} = curBlock(:,4);
%             xy2_data2{blockIdx,s} = curBlock(:,5);
            idx = 1;
            for x = 1:length(xy1_data1{blockIdx,s-3})-1
%             for x = 1:length(xy1_data1{blockIdx,s})-1
                if abs(xy1_data1{blockIdx,s-3}(x+1) - xy1_data1{blockIdx,s-3}(x)) > 200
%                 if abs(xy1_data1{blockIdx,s}(x+1) - xy1_data1{blockIdx,s}(x)) > 200
                    trialIdx(idx) = x; % gets the last sample of the trial
                    idx = idx + 1;
                end
            end
            trialIdx(idx) = length(xy1_data1{blockIdx,s-3});
            trialEndpoint(:,blockIdx,s-3) = trialIdx;
%             trialIdx(idx) = length(xy1_data1{blockIdx,s});
%             trialEndpoint(:,blockIdx,s) = trialIdx;
            blockIdx = blockIdx + 1;
            
%             final_xy1_d1(:,blockIdx,s) = xy1_data1{blockIdx,s}(trialIdx);
%             final_xy1_d2(:,blockIdx,s) = xy1_data2{blockIdx,s}(trialIdx);
%             final_xy2_d1(:,blockIdx,s) = xy2_data1{blockIdx,s}(trialIdx);
%             final_xy2_d2(:,blockIdx,s) = xy2_data2{blockIdx,s}(trialIdx);
%             
%             RH_xy1_d1(:,blockIdx,s) = final_xy1_d1(RH_idx,blockIdx,s); % len = 48 bc 48 RH trials per 64-trial block
%             RH_xy1_d2(:,blockIdx,s) = final_xy1_d2(RH_idx,blockIdx,s);
%             RH_xy2_d1(:,blockIdx,s) = final_xy2_d1(RH_idx,blockIdx,s);
%             RH_xy2_d2(:,blockIdx,s) = final_xy2_d2(RH_idx,blockIdx,s);
        end
% 
%         % RT ANALYSIS (for all trials)
% %         S1_onset_time = data.S1_onset_time;
% %         S2_onset_time = data.S2_onset_time;
        rt = data.rt;

        if ~data.grasping
           nInvalidKey(s,idx) = length(find(rt<(.1 + .3 + .1 + .1))); 
           idx = idx + 1;
        else
           acc_grasp = data.acc_grasp;
           SOT_data = data.SOT_data;
           percAccGrasp(s,grIdx) = sum(data.acc_grasp)/length(data.acc_grasp);
%            initLatency = SOT_data(end) - S1_onset_time 
           grIdx = grIdx + 1;
        end
        
    end % END OF BLOCKWISE ANALYSIS
    
    if s >= 4
        angleFinal(:,b-1) = data.angleGraspFinal;
        for ii = 1:8
            NC_angle = angleFinal(NC_idx{ii,s},:);
            C_angle = angleFinal(C_idx{ii,s},:);
            RH_NC_angle = angleFinal([RH_H_NC_idx{ii,s} RH_V_NC_idx{ii,s}],:);
            RH_C_angle = angleFinal([RH_H_C_idx{ii,s} RH_V_C_idx{ii,s}],:);
            RH_V_NC_angle = angleFinal(RH_V_NC_idx{ii,s},:);
            RH_V_C_angle = angleFinal(RH_V_C_idx{ii,s},:);
        end
        NC_angle = [NC_angle(:,1);NC_angle(:,3);NC_angle(:,5);NC_angle(:,7)];
        C_angle = [C_angle(:,1);C_angle(:,3);C_angle(:,5);C_angle(:,7)];
        RH_NC_angle = [RH_NC_angle(:,1);RH_NC_angle(:,3);RH_NC_angle(:,5);RH_NC_angle(:,7)];
        RH_C_angle = [RH_C_angle(:,1);RH_C_angle(:,3);RH_C_angle(:,5);RH_C_angle(:,7)];
        RH_V_NC_angle = [RH_V_NC_angle(:,1);RH_V_NC_angle(:,3);RH_V_NC_angle(:,5);RH_V_NC_angle(:,7)];
        RH_V_C_angle = [RH_V_C_angle(:,1);RH_V_C_angle(:,3);RH_V_C_angle(:,5);RH_V_C_angle(:,7)];
    end

    RH_diff(s) = length(find(~RH_NC_Ans_key)); % = 48, because 48 total DIFF trials per task type (grasp/key) per subj; 12 * 4
    RH_same(s) = length(find(RH_NC_Ans_key));
    LH_diff(s) = length(find(~LH_NC_Ans_key));% = 16, because 16 total DIFF trials per task type per subj; 4 * 4
    LH_same(s) = length(find(LH_NC_Ans_key));
    
    [RH_NC_data_key(s), RH_NC_dp_key(s)] = dprimeCalc(RH_NC_Ans_key,RH_NC_Rsp_key,RH_same(s),RH_diff(s));
    [RH_NC_data_grp(s), RH_NC_dp_grp(s)] = dprimeCalc(RH_NC_Ans_grp,RH_NC_Rsp_grp,RH_same(s),RH_diff(s));
    [RH_C_data_key(s), RH_C_dp_key(s)] = dprimeCalc(RH_C_Ans_key,RH_C_Rsp_key,RH_same(s),RH_diff(s));
    [RH_C_data_grp(s), RH_C_dp_grp(s)] = dprimeCalc(RH_C_Ans_grp,RH_C_Rsp_grp,RH_same(s),RH_diff(s));
    
    [LH_NC_data_key(s), LH_NC_dp_key(s)] = dprimeCalc(LH_NC_Ans_key,LH_NC_Rsp_key,LH_same(s),LH_diff(s));
    [LH_NC_data_grp(s), LH_NC_dp_grp(s)] = dprimeCalc(LH_NC_Ans_grp,LH_NC_Rsp_grp,LH_same(s),LH_diff(s));
    [LH_C_data_key(s), LH_C_dp_key(s)] = dprimeCalc(LH_C_Ans_key,LH_C_Rsp_key,LH_same(s),LH_diff(s));
    [LH_C_data_grp(s), LH_C_dp_grp(s)] = dprimeCalc(LH_C_Ans_grp,LH_C_Rsp_grp,LH_same(s),LH_diff(s));
end

mean_RH_NC_key = mean(RH_NC_dp_key);
mean_RH_C_key = mean(RH_C_dp_key);
mean_RH_NC_grp = mean(RH_NC_dp_grp);
mean_RH_C_grp = mean(RH_C_dp_grp);

err_RH_NC_key = std(RH_NC_dp_key)/sqrt(length(RH_NC_dp_key));
err_RH_C_key = std(RH_C_dp_key)/sqrt(length(RH_NC_dp_key));
err_RH_NC_grp = std(RH_NC_dp_grp)/sqrt(length(RH_NC_dp_key));
err_RH_C_grp = std(RH_C_dp_grp)/sqrt(length(RH_NC_dp_key));

mean_RH_NC_C_key = mean(RH_NC_dp_key-RH_C_dp_key);
mean_RH_NC_C_grp = mean(RH_NC_dp_grp-RH_C_dp_grp);
err_RH_NC_C_key = std(RH_NC_dp_key-RH_C_dp_key)/sqrt(length(RH_NC_dp_key-RH_C_dp_key));
err_RH_NC_C_grp = std(RH_NC_dp_grp-RH_C_dp_grp)/sqrt(length(RH_NC_dp_key-RH_C_dp_key));

% COMPARE GRASP/KEY D-PRIMES
RH_NC_diff = RH_NC_dp_grp - RH_NC_dp_key;
RH_C_diff = RH_C_dp_grp - RH_C_dp_key;
LH_NC_diff = LH_NC_dp_grp - LH_NC_dp_key;
LH_C_diff = LH_C_dp_grp - LH_C_dp_key;

%% SENSITIVITY PLOT
cndNames = {'No crowding','Crowding'};
x1 = [1 1.75];
% gcaOpts = {'XTick',x1,'XTickLabel',cndNames,'box','off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};
gcaOpts = {'XTick',[],'box','off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};
figure

x = [1./mean(thresholds(:,1)) 0; 1./mean(thresholds(:,2)) 1];
err = [std(1./thresholds(:,1))/sqrt(length(thresholds(:,1))); std(1./thresholds(:,2))/sqrt(length(thresholds(:,2)))];

for i=1:length(x)
    h = bar(x1(i),x(i,1),0.5,'EdgeColor',[0 0 0],'LineWidth',1.5);
    if i == 1, hold on, end
    if ~x(i,2)
        col = [12,121,153]/255;
    else
        col = [240,59,37]/255;
    end
    set(h, 'FaceColor', col, 'FaceAlpha', 0.5)
end
hold on;
ngroups = 1;
nbars = size(x, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    errorbar(x1(i), x(i), err(i), 'k','LineWidth',1.5,'linestyle', 'none');
end
legend(cndNames,'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Orientation discrimination sensitivity')
xlabel('Condition')
ylabel('Sensitivity')
ylim([0 .7])

%% PLOT D PRIMES / WITH GRASP / FOR RIGHT HEMI
cndNames = {'No crowding','Crowding'};
xtick = [1.375 3.375];
x1 = [1 1.75 3 3.75];
gcaOpts = {'XTick',xtick,'XTickLabels',cndNames,'box','off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
xdata = [mean_RH_NC_key, 1, 0; 
    mean_RH_NC_grp, 2, 1; 
    mean_RH_C_key, 3, 2; 
    mean_RH_C_grp, 4, 3];
err = [std(RH_NC_dp_key)/sqrt(nSubj); std(RH_NC_dp_grp)/sqrt(nSubj); std(RH_C_dp_key)/sqrt(nSubj); std(RH_C_dp_grp)/sqrt(nSubj)];

for i=1:length(xdata)
    h = bar(x1(i),xdata(i,1),0.7,'EdgeColor',[0 0 0],'LineWidth',1.5);
    if i == 1, hold on, end
    if ~xdata(i,3)
        col = [12,121,153]/255;
        transparency = 0.5;
    elseif xdata(i,3) == 1
        col = [12,121,153]/255; 
        transparency = 1;
    elseif xdata(i,3) == 2
        col = [240,59,37]/255;
        transparency = 0.5;
    else
        col = [240,59,37]/255;
        transparency = 1;
    end
    set(h, 'FaceColor', col, 'FaceAlpha',transparency)
end
g(1) = bar(NaN,NaN,'FaceColor', 'k', 'FaceAlpha',0.25);
g(2) = bar(NaN,NaN,'FaceColor', 'k', 'FaceAlpha',.75);
hold on;
ngroups = 1;
nbars = size(xdata, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    errorbar(x1(i), xdata(i), err(i), 'k', 'LineWidth',1.5,'linestyle', 'none');
end
legend(g,{'Keypress-only','Grasp'},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('Orientation sensitivity (d'')')
xlabel('Condition')
ylabel('d-prime')
ylim([0 2.7])

% figure
% h = bar([RH_NC_dp_key, mean_RH_NC_key; RH_NC_dp_grp, mean_RH_NC_grp; RH_C_dp_key, mean_RH_C_key; RH_C_dp_grp, mean_RH_C_grp]); 
% h(1).FaceColor = [135,205,215]/255;
% h(2).FaceColor = [240,59,37]/255;
% h(3).FaceColor = [250,185,219]/255;
% h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% % h(6).FaceColor = [0 0 0]/255;
% legend({subj{:},'avg'},'AutoUpdate','off','location','northeast');
% 
% set(gca,gcaOpts{:})
% title('Orientation sensitivity (d'')')
% xlabel('Condition')
% ylabel('d-prime')
% ylim([-1 3.5])

%% PLOT D PRIME COMPARISON, GRASP VS NO GRASP / RIGHT HEMI
% cndNames = {'No crowding','Crowding'};
% gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
%     'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};
% 
% figure
% % h = bar([RH_NC_diff, mean_RH_NC_grp-mean_RH_NC_key; RH_C_diff, mean_RH_C_grp-mean_RH_C_key]); 
% % h(1).FaceColor = [135,205,215]/255;
% % h(2).FaceColor = [240,59,37]/255;
% % h(3).FaceColor = [250,185,219]/255;
% % h(4).FaceColor = [28 15 142]/255;
% % h(5).FaceColor = [255 255 255]/255;
% h = bar([mean_RH_NC_grp-mean_RH_NC_key; mean_RH_C_grp-mean_RH_C_key]); 
% h(1).FaceColor = [135,205,215]/255;
% % h(6).FaceColor = [0 0 0]/255;
% legend({'avg'},'AutoUpdate','off','location','northeast');
% 
% set(gca,gcaOpts{:})
% title('d'' diff (grasping - key only) in RH')
% xlabel('Condition')
% ylabel('d'' difference')
% ylim([-1 1])

cndNames = {'NC','C'};
x1 = [1 1.75];
% gcaOpts = {'XTick',x1,'XTickLabel',cndNames,'box','off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};
gcaOpts = {'XTick',[],'box','off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
x = [mean(RH_NC_dp_grp-RH_NC_dp_key) 0; mean(RH_C_dp_grp-RH_C_dp_key) 1];
err = [std(RH_NC_dp_grp-RH_NC_dp_key)/sqrt(nSubj); std(RH_C_dp_grp-RH_C_dp_key)/sqrt(nSubj)];

for i=1:length(x)
    h = bar(x1(i),x(i,1),0.5,'EdgeColor',[0 0 0],'LineWidth',1.5);
    if i == 1, hold on, end
    if ~x(i,2)
        col = [12,121,153]/255;
    elseif x(i,2) == 1
        col = [240,59,37]/255;
    end
    set(h, 'FaceColor', col, 'FaceAlpha', 0.75)
end
hold on;
ngroups = 1;
nbars = size(x, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    errorbar(x1(i), x(i), err(i), 'k', 'LineWidth', 1.5, 'linestyle', 'none');
end
line([0 2.75],[0 0],'LineWidth',1.5,'color','k')
legend(cndNames,'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Grasp-keypress sensitivity difference')
xlabel('Condition')
ylabel('d'' difference')
xlim([.25 2.5])
ylim([-.4 1])

%% PLOT GRASP ANGLE ERROR

avg_NC_angle = mean(NC_angle);
avg_C_angle = mean(C_angle);
avg_RH_NC_angle = mean(RH_NC_angle);
avg_RH_C_angle = mean(RH_C_angle);
avg_RH_V_NC_angle = mean(RH_V_NC_angle);
avg_RH_V_C_angle = mean(RH_V_C_angle);

err_NC_angle = std(NC_angle)/sqrt(length(NC_angle));
err_C_angle = std(C_angle)/sqrt(length(C_angle));
err_RH_NC_angle = std(RH_NC_angle)/sqrt(length(RH_NC_angle));
err_RH_C_angle = std(RH_C_angle)/sqrt(length(RH_C_angle));
err_RH_V_NC_angle = std(RH_V_NC_angle)/sqrt(length(RH_V_NC_angle));
err_RH_V_C_angle = std(RH_V_C_angle)/sqrt(length(RH_V_C_angle));
    
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
% x = [avg_NC_angle; avg_C_angle];
% err = [err_NC_
% x = [avg_RH_NC_angle; avg_RH_C_angle];
% err = [err_RH_NC_angle; err_RH_C_angle];
x = [avg_RH_V_NC_angle; avg_RH_V_C_angle];
err = [err_RH_V_NC_angle; err_RH_V_C_angle];

h = bar(x);
h(1).FaceColor = [135,205,215]/255;
hold on;
ngroups = 1;
nbars = size(x, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x1 = i;
    errorbar(x1, x(i), err(i), 'k', 'linestyle', 'none');
end
legend({'Average','SEM'},'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Grasp angle error')
xlabel('Condition')
ylabel('Error (degrees)')
ylim([0 9])

%% PLOT THRESHOLDS
close all
cndNames = {'No crowding','Crowding'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
x = [mean(thresholds(:,1)); mean(thresholds(:,2))];
err = [std(thresholds(:,1))/sqrt(length(thresholds(:,1)));std(thresholds(:,2))/sqrt(length(thresholds(:,2)))];

h = bar(x);
h(1).FaceColor = [135,205,215]/255;
hold on;
ngroups = 1;
nbars = size(x, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x1 = i;
%     x1 = (1:ngroups) - (groupwidth/2) + (2*i-1) * (groupwidth / (2*nbars));
    errorbar(x1, x(i), err(i), 'k', 'linestyle', 'none');
end
legend({'Average','SEM'},'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Orientation discrimination thresholds, 80% accuracy')
% xlabel('Condition')
ylabel('80% threshold (degrees)')
ylim([0 6])


%%

figure
x1 = thresholds(:,1)';
x2 = thresholds(:,2)';
h = bar([x1;x2]);
h(1).FaceColor = [135,205,215]/255;
h(2).FaceColor = [240,59,37]/255;
h(3).FaceColor = [250,185,219]/255;
h(4).FaceColor = [28 15 142]/255;
% h(5).FaceColor = [255 255 255]/255;
% h(6).FaceColor = [0 0 0]/255;

legend({subj{:}},'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Orientation thresholds')
xlabel('Crowding condition')
ylabel('80% accuracy threshold (degrees)')
ylim([0 6])

%% PLOT TRAJECTORIES
horiz_xy1_d1 = [];horiz_xy1_d2 = [];horiz_xy2_d1 = [];horiz_xy2_d2 = [];
vert_xy1_d1 = [];vert_xy1_d2 = [];vert_xy2_d1 = [];vert_xy2_d2 = [];
H_NC_xy1_d1 = [];H_NC_xy1_d2 = [];H_NC_xy2_d1 = [];H_NC_xy2_d2 = [];
V_NC_xy1_d1 = [];V_NC_xy1_d2 = [];V_NC_xy2_d1 = [];V_NC_xy2_d2 = [];
H_C_xy1_d1 = [];H_C_xy1_d2 = [];H_C_xy2_d1 = [];H_C_xy2_d2 = [];
V_C_xy1_d1 = [];V_C_xy1_d2 = [];V_C_xy2_d1 = [];V_C_xy2_d2 = [];
 
% grp_horizIdx = horizIdx(1:32,find(grasping(:,4)),4);
% grp_vertIdx = vertIdx(1:32,find(grasping(:,4)),4);

% grp_RH_H_NC_idx = RH_H_NC_idx(find(grasping(:,4)),4); % eventually replace '4' with subj loop
% grp_H_NC_idx = [RH_H_NC_idx(find(grasping(:,4)),4);LH_H_NC_idx(find(grasping(:,4)),4)];

grp_H_NC_idx = H_NC_idx(find(grasping(:,s)),s);
grp_H_C_idx = H_C_idx(find(grasping(:,s)),s);
grp_V_NC_idx = V_NC_idx(find(grasping(:,s)),s);
grp_V_C_idx = V_C_idx(find(grasping(:,s)),s);

nBlock = size(xy1_data1,1);

split_xy1_d1 = cell(allTrials(end),nBlock,nSubj-3);
split_xy1_d2 = cell(allTrials(end),nBlock,nSubj-3);
split_xy2_d1 = cell(allTrials(end),nBlock,nSubj-3);
split_xy2_d2 = cell(allTrials(end),nBlock,nSubj-3);

% split trajectories into trials
for s = 1:nSubj-3
    for b = 1:size(xy1_data1,1)
        split_xy1_d1{1,b,s} = xy1_data1{b,s}(1:trialEndpoint(1,b,s));
        split_xy1_d2{1,b,s} = xy1_data2{b,s}(1:trialEndpoint(1,b,s))*-1;
        split_xy2_d1{1,b,s} = xy2_data1{b,s}(1:trialEndpoint(1,b,s));
        split_xy2_d2{1,b,s} = xy2_data2{b,s}(1:trialEndpoint(1,b,s))*-1;
        for u = 2:length(trialEndpoint(:,b,s))
            split_xy1_d1{u,b,s} = xy1_data1{b,s}(trialEndpoint(u-1,b,s)+1:trialEndpoint(u,b,s));
            split_xy1_d2{u,b,s} = xy1_data2{b,s}(trialEndpoint(u-1,b,s)+1:trialEndpoint(u,b,s))*-1;
            split_xy2_d1{u,b,s} = xy2_data1{b,s}(trialEndpoint(u-1,b,s)+1:trialEndpoint(u,b,s));
            split_xy2_d2{u,b,s} = xy2_data2{b,s}(trialEndpoint(u-1,b,s)+1:trialEndpoint(u,b,s))*-1;
        end
    end
end

for i = 1:size(split_xy1_d1,2)    
    H_NC_xy1_d1 = [H_NC_xy1_d1; split_xy1_d1(grp_H_NC_idx{i},i)];
    H_NC_xy1_d2 = [H_NC_xy1_d2; split_xy1_d2(grp_H_NC_idx{i},i)];
    H_NC_xy2_d1 = [H_NC_xy2_d1; split_xy2_d1(grp_H_NC_idx{i},i)];
    H_NC_xy2_d2 = [H_NC_xy2_d2; split_xy2_d2(grp_H_NC_idx{i},i)];
    
    V_NC_xy1_d1 = [V_NC_xy1_d1; split_xy1_d1(grp_V_NC_idx{i},i)];
    V_NC_xy1_d2 = [V_NC_xy1_d2; split_xy1_d2(grp_V_NC_idx{i},i)];
    V_NC_xy2_d1 = [V_NC_xy2_d1; split_xy2_d1(grp_V_NC_idx{i},i)];
    V_NC_xy2_d2 = [V_NC_xy2_d2; split_xy2_d2(grp_V_NC_idx{i},i)];
    
    H_C_xy1_d1 = [H_C_xy1_d1; split_xy1_d1(grp_H_C_idx{i},i)];
    H_C_xy1_d2 = [H_C_xy1_d2; split_xy1_d2(grp_H_C_idx{i},i)];
    H_C_xy2_d1 = [H_C_xy2_d1; split_xy2_d1(grp_H_C_idx{i},i)];
    H_C_xy2_d2 = [H_C_xy2_d2; split_xy2_d2(grp_H_C_idx{i},i)];
    
    V_C_xy1_d1 = [V_C_xy1_d1; split_xy1_d1(grp_V_C_idx{i},i)];
    V_C_xy1_d2 = [V_C_xy1_d2; split_xy1_d2(grp_V_C_idx{i},i)];
    V_C_xy2_d1 = [V_C_xy2_d1; split_xy2_d1(grp_V_C_idx{i},i)];
    V_C_xy2_d2 = [V_C_xy2_d2; split_xy2_d2(grp_V_C_idx{i},i)];
end

%% PLOTTING TRAJECTORIES  
gcaOpts = {'box','off','XTick',[],'YTick',[],'tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

% figure
% for s = 1:nSubj-3
%     subplot(2,nSubj-3,s)
%     hold on
%     for i = 1:length(horiz_xy1_d1)
%         plot(horiz_xy1_d1{i},horiz_xy1_d2{i},'color','r')
%         plot(horiz_xy2_d1{i},horiz_xy2_d2{i},'color','b')
%     end
%     set(gca,gcaOpts{:})
%     ylim([-1250 -300]);
%     xlim([-200 1400]);
%     title('All horizontal trajectories')
%     
%     subplot(2,nSubj-3,s+nSubj-3)
%     hold on
%     for i = 1:length(vert_xy1_d1)
%         plot(vert_xy1_d1{i},vert_xy1_d2{i},'color','r')
%         plot(vert_xy2_d1{i},vert_xy2_d2{i},'color','b')
%     end
%     set(gca,gcaOpts{:})
%     ylim([-1250 -300]);
%     xlim([-200 1400]);
%     title('All vertical trajectories')
% end

figure
for s = 1:nSubj-3
    subplot(2,2,1)
    hold on
    for i = 1:length(H_NC_xy1_d1)
        plot(H_NC_xy1_d1{i},H_NC_xy1_d2{i},'color','r')
        plot(H_NC_xy2_d1{i},H_NC_xy2_d2{i},'color','b')
        xy1d1(i) = H_NC_xy1_d1{i}(end);
        xy1d2(i) = H_NC_xy1_d2{i}(end);
        xy2d1(i) = H_NC_xy2_d1{i}(end);
        xy2d2(i) = H_NC_xy2_d2{i}(end);
    end
    RH = find(xy1d1 > 400);
    LH = find(xy1d1 <= 400);
    H_NC_end_xy1d1 = mean(xy1d1(RH));
    H_NC_end_xy1d2 = mean(xy1d2(RH));
    H_NC_end_xy2d1 = mean(xy2d1(RH));
    H_NC_end_xy2d2 = mean(xy2d2(RH));
    H_NC_err_xy1d1 = std(xy1d1(RH))/sqrt(length(xy1d1(RH)));
    H_NC_err_xy1d2 = std(xy1d2(RH))/sqrt(length(xy1d2(RH)));
    H_NC_err_xy2d1 = std(xy2d1(RH))/sqrt(length(xy2d1(RH)));
    H_NC_err_xy2d2 = std(xy2d2(RH))/sqrt(length(xy2d2(RH)));
    LH_H_NC_end_xy1d1 = mean(xy1d1(LH));
    LH_H_NC_end_xy1d2 = mean(xy1d2(LH));
    LH_H_NC_end_xy2d1 = mean(xy2d1(LH));
    LH_H_NC_end_xy2d2 = mean(xy2d2(LH));
    LH_H_NC_err_xy1d1 = std(xy1d1(LH))/sqrt(length(xy1d1(LH)));
    LH_H_NC_err_xy1d2 = std(xy1d2(LH))/sqrt(length(xy1d2(LH)));
    LH_H_NC_err_xy2d1 = std(xy2d1(LH))/sqrt(length(xy2d1(LH)));
    LH_H_NC_err_xy2d2 = std(xy2d2(LH))/sqrt(length(xy2d2(LH)));
    clear xy1d1; clear xy1d2; clear xy2d1; clear xy2d2;
    
    plot([H_NC_end_xy1d1 H_NC_end_xy2d1],[H_NC_end_xy1d2 H_NC_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8)
    plot([LH_H_NC_end_xy1d1 LH_H_NC_end_xy2d1],[LH_H_NC_end_xy1d2 LH_H_NC_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8)
    set(gca,gcaOpts{:})
    ylim([-1250 -300]);
    xlim([-200 1400]);
    pbaspect([1 1 1])
    title('Non-crowded horiz. grasp')
    
    subplot(2,2,2)
    hold on
    for i = 1:length(V_NC_xy1_d1)
        plot(V_NC_xy1_d1{i},V_NC_xy1_d2{i},'color','r')
        plot(V_NC_xy2_d1{i},V_NC_xy2_d2{i},'color','b')
        xy1d1(i) = V_NC_xy1_d1{i}(end);
        xy1d2(i) = V_NC_xy1_d2{i}(end);
        xy2d1(i) = V_NC_xy2_d1{i}(end);
        xy2d2(i) = V_NC_xy2_d2{i}(end);
    end
    RH = find(xy1d1 > 400);
    LH = find(xy1d1 <= 400);
    V_NC_end_xy1d1 = mean(xy1d1(RH));
    V_NC_end_xy1d2 = mean(xy1d2(RH));
    V_NC_end_xy2d1 = mean(xy2d1(RH));
    V_NC_end_xy2d2 = mean(xy2d2(RH));
    V_NC_err_xy1d1 = std(xy1d1(RH))/sqrt(length(xy1d1(RH))); % calculate variance 
    V_NC_err_xy1d2 = std(xy1d2(RH))/sqrt(length(xy1d2(RH)));
    V_NC_err_xy2d1 = std(xy2d1(RH))/sqrt(length(xy1d2(RH)));
    V_NC_err_xy2d2 = std(xy2d2(RH))/sqrt(length(xy2d2(RH)));
    LH_V_NC_end_xy1d1 = mean(xy1d1(LH));
    LH_V_NC_end_xy1d2 = mean(xy1d2(LH));
    LH_V_NC_end_xy2d1 = mean(xy2d1(LH));
    LH_V_NC_end_xy2d2 = mean(xy2d2(LH));
    LH_V_NC_err_xy1d1 = std(xy1d1(LH))/sqrt(length(xy1d1(LH)));
    LH_V_NC_err_xy1d2 = std(xy1d2(LH))/sqrt(length(xy1d2(LH)));
    LH_V_NC_err_xy2d1 = std(xy2d1(LH))/sqrt(length(xy2d1(LH)));
    LH_V_NC_err_xy2d2 = std(xy2d2(LH))/sqrt(length(xy2d2(LH)));
    clear xy1d1; clear xy1d2; clear xy2d1; clear xy2d2;
    
    plot([V_NC_end_xy1d1 V_NC_end_xy2d1],[V_NC_end_xy1d2 V_NC_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8)
    plot([LH_V_NC_end_xy1d1 LH_V_NC_end_xy2d1],[LH_V_NC_end_xy1d2 LH_V_NC_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8)
    set(gca,gcaOpts{:})
    ylim([-1250 -300]);
    xlim([-200 1400]);
    pbaspect([1 1 1])
    title('Non-crowded vert. grasp')
end


for s = 1:nSubj-3
    subplot(2,2,3)
    hold on
    for i = 1:length(H_C_xy1_d1)
        plot(H_C_xy1_d1{i},H_C_xy1_d2{i},'color','r')
        plot(H_C_xy2_d1{i},H_C_xy2_d2{i},'color','b')
        xy1d1(i) = H_C_xy1_d1{i}(end);
        xy1d2(i) = H_C_xy1_d2{i}(end);
        xy2d1(i) = H_C_xy2_d1{i}(end);
        xy2d2(i) = H_C_xy2_d2{i}(end);
    end
    RH = find(xy1d1 > 400);
    LH = find(xy1d1 <= 400);
    H_C_end_xy1d1 = mean(xy1d1(RH));
    H_C_end_xy1d2 = mean(xy1d2(RH));
    H_C_end_xy2d1 = mean(xy2d1(RH));
    H_C_end_xy2d2 = mean(xy2d2(RH));
    H_C_err_xy1d1 = std(xy1d1(RH))/sqrt(length(xy1d1(RH)));
    H_C_err_xy1d2 = std(xy1d2(RH))/sqrt(length(xy1d2(RH)));
    H_C_err_xy2d1 = std(xy2d1(RH))/sqrt(length(xy1d2(RH)));
    H_C_err_xy2d2 = std(xy2d2(RH))/sqrt(length(xy2d2(RH)));
    LH_H_C_end_xy1d1 = mean(xy1d1(LH));
    LH_H_C_end_xy1d2 = mean(xy1d2(LH));
    LH_H_C_end_xy2d1 = mean(xy2d1(LH));
    LH_H_C_end_xy2d2 = mean(xy2d2(LH));
    LH_H_C_err_xy1d1 = std(xy1d1(LH))/sqrt(length(xy1d1(LH)));
    LH_H_C_err_xy1d2 = std(xy1d2(LH))/sqrt(length(xy1d2(LH)));
    LH_H_C_err_xy2d1 = std(xy2d1(LH))/sqrt(length(xy1d2(LH)));
    LH_H_C_err_xy2d2 = std(xy2d2(LH))/sqrt(length(xy2d2(LH)));
    clear xy1d1; clear xy1d2; clear xy2d1; clear xy2d2;
    
    plot([H_C_end_xy1d1 H_C_end_xy2d1],[H_C_end_xy1d2 H_C_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8)
    plot([LH_H_C_end_xy1d1 LH_H_C_end_xy2d1],[LH_H_C_end_xy1d2 LH_H_C_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8)
    set(gca,gcaOpts{:})
    ylim([-1250 -300]);
    xlim([-200 1400]);
    pbaspect([1 1 1])
    title('Crowded horiz. grasp')
    
    subplot(2,2,4)
    hold on
    for i = 1:length(V_C_xy1_d1)
        H(i,1) = plot(V_C_xy1_d1{i},V_C_xy1_d2{i},'color','r');
        H(i,2) = plot(V_C_xy2_d1{i},V_C_xy2_d2{i},'color','b');
        xy1d1(i) = V_C_xy1_d1{i}(end);
        xy1d2(i) = V_C_xy1_d2{i}(end);
        xy2d1(i) = V_C_xy2_d1{i}(end);
        xy2d2(i) = V_C_xy2_d2{i}(end);
    end
    RH = find(xy1d1 > 400);
    LH = find(xy1d1 <= 400);
    V_C_end_xy1d1 = mean(xy1d1(RH));
    V_C_end_xy1d2 = mean(xy1d2(RH));
    V_C_end_xy2d1 = mean(xy2d1(RH));
    V_C_end_xy2d2 = mean(xy2d2(RH));
    V_C_err_xy1d1 = std(xy1d1(RH))/sqrt(length(xy1d1(RH)));
    V_C_err_xy1d2 = std(xy1d2(RH))/sqrt(length(xy1d2(RH)));
    V_C_err_xy2d1 = std(xy2d1(RH))/sqrt(length(xy1d2(RH)));
    V_C_err_xy2d2 = std(xy2d2(RH))/sqrt(length(xy2d2(RH)));
    LH_V_C_end_xy1d1 = mean(xy1d1(LH));
    LH_V_C_end_xy1d2 = mean(xy1d2(LH));
    LH_V_C_end_xy2d1 = mean(xy2d1(LH));
    LH_V_C_end_xy2d2 = mean(xy2d2(LH));
    LH_V_C_err_xy1d1 = std(xy1d1(LH))/sqrt(length(xy1d1(LH)));
    LH_V_C_err_xy1d2 = std(xy1d2(LH))/sqrt(length(xy1d2(LH)));
    LH_V_C_err_xy2d1 = std(xy2d1(LH))/sqrt(length(xy1d2(LH)));
    LH_V_C_err_xy2d2 = std(xy2d2(LH))/sqrt(length(xy2d2(LH)));
    clear xy1d1; clear xy1d2; clear xy2d1; clear xy2d2;
    
    plot([V_C_end_xy1d1 V_C_end_xy2d1],[V_C_end_xy1d2 V_C_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8);
    plot([LH_V_C_end_xy1d1 LH_V_C_end_xy2d1],[LH_V_C_end_xy1d2 LH_V_C_end_xy2d2],'-o','color','c','LineWidth',2,'MarkerSize',8)
    legend({'index','thumb'},'AutoUpdate','off','location','southeast');
    set(gca,gcaOpts{:})
    xlabel('x position')
    ylabel('y position')
    ylim([-1250 -300]);
    xlim([-200 1400]);
    pbaspect([1 1 1])
    title('Crowded vert. grasp')
end

% H_NC = [H_NC_err_xy1d1 H_NC_err_xy1d2]; %  H_NC_err_xy2d1 H_NC_err_xy2d2
% V_NC = [V_NC_err_xy1d1 V_NC_err_xy1d2]; %  V_NC_err_xy2d1 V_NC_err_xy2d2
% H_C = [H_C_err_xy1d1 H_C_err_xy1d2]; %  H_C_err_xy2d1 H_C_err_xy2d2
% V_C = [V_C_err_xy1d1 V_C_err_xy1d2]; %  V_C_err_xy2d1 V_C_err_xy2d2
% 
% gcaOpts = {'XTick',1:4,'XTickLabel',{'Horiz, NC' 'Vert, NC' 'Horiz, C' 'Vert, C'},'box',...
%     'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};
% figure
% bar([mean(H_NC) mean(V_NC) mean(H_C) mean(V_C)])
% set(gca,gcaOpts{:})
% title('Standard errors for each endpoint')
% % legend({'Horiz, NC' 'Vert, NC' 'Horiz, C' 'Vert, C'})


%% FILTERED TRAJECTORY PLOTS

% C_trialEndpoint = trialEndpoint(C_idx,:,:);
% figure
% for s = 1:nSubj
%     subplot(1,nSubj,s)
%     hold on
%     for b = 1:4
%         plot(xy1_data1{b,s}(1:C_trialEndpoint(1,b,s)),xy1_data2{b,s}(1:C_trialEndpoint(1,b,s))*-1,'b')
%         plot(xy2_data1{b,s}(1:C_trialEndpoint(1,b,s)),xy2_data2{b,s}(1:C_trialEndpoint(1,b,s))*-1,'r')
%         for u = 1:length(C_trialEndpoint(:,b,s))-1
%             plot(xy1_data1{b,s}(C_trialEndpoint(u,b,s)+1:C_trialEndpoint(u+1,b,s)),xy1_data2{b,s}(C_trialEndpoint(u,b,s)+...
%                 1:C_trialEndpoint(u+1,b,s))*-1,'color','b')
%             plot(xy2_data1{b,s}(C_trialEndpoint(u,b,s)+1:C_trialEndpoint(u+1,b,s)),xy2_data2{b,s}(C_trialEndpoint(u,b,s)+...
%                 1:C_trialEndpoint(u+1,b,s))*-1,'color','r')
%         end
%     end
% end
% set(gca,gcaOpts{:})
% ylim([-1250 -300]);
% xlim([-200 1400]);
% title('Trajectories for crowding trials')
% 
% NC_trialEndpoint = trialEndpoint(NC_idx,:,:);
% figure
% for s = 1:nSubj
%     subplot(1,nSubj,s)
%     hold on
%     for b = 1:4
%         plot(xy1_data1{b,s}(1:NC_trialEndpoint(1,b,s)),xy1_data2{b,s}(1:NC_trialEndpoint(1,b,s))*-1,'b')
%         plot(xy2_data1{b,s}(1:NC_trialEndpoint(1,b,s)),xy2_data2{b,s}(1:NC_trialEndpoint(1,b,s))*-1,'r')
%         for u = 1:length(NC_trialEndpoint(:,b,s))-1
%             plot(xy1_data1{b,s}(NC_trialEndpoint(u,b,s)+1:NC_trialEndpoint(u+1,b,s)),xy1_data2{b,s}(NC_trialEndpoint(u,b,s)+...
%                 1:NC_trialEndpoint(u+1,b,s))*-1,'color','b')
%             plot(xy2_data1{b,s}(NC_trialEndpoint(u,b,s)+1:NC_trialEndpoint(u+1,b,s)),xy2_data2{b,s}(NC_trialEndpoint(u,b,s)+...
%                 1:NC_trialEndpoint(u+1,b,s))*-1,'color','r')
%         end
%     end
% end
% set(gca,gcaOpts{:})
% ylim([-1250 -300]);
% xlim([-200 1400]);
% title('Trajectories for non-crowding trials')

%% PLOT HITS etc FOR RIGHT HEMI
cndNames = {'H,K','M,K','F,K','C,K','H,G','M,G','F,G','C,G'};
gcaOpts = {'XTick',1:8,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',9};

figure
for i = 1:nSubj
    subplot(2,nSubj,i)
    bar([RH_NC_data_key(i).hit; RH_NC_data_key(i).miss; RH_NC_data_key(i).falseAlarm; RH_NC_data_key(i).corReject; ...
    RH_NC_data_grp(i).hit; RH_NC_data_grp(i).miss; RH_NC_data_grp(i).falseAlarm; RH_NC_data_grp(i).corReject]) 
    set(gca,gcaOpts{:})
    ylim([0 50])
    title(sprintf('RH non-crowd S0%d',i))
    xlabel('Condition')
    ylabel('d-prime')
    
    subplot(2,nSubj,i+nSubj)
    bar([RH_C_data_key(i).hit; RH_C_data_key(i).miss; RH_C_data_key(i).falseAlarm; RH_C_data_key(i).corReject; ...
    RH_C_data_grp(i).hit; RH_C_data_grp(i).miss; RH_C_data_grp(i).falseAlarm; RH_C_data_grp(i).corReject]) 
    set(gca,gcaOpts{:})
    ylim([0 50])
    title(sprintf('RH crowd S0%d',i))
    xlabel('Condition')
    ylabel('d-prime')
end



%% COMPARE C VS NC DPRIMES 
cndNames = {'Grasp','Key-only'};
gcaOpts = {'XTick',[],'box','off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

% figure
x1 = [1 1.75];
xdata = [mean_RH_NC_C_grp, 1, 1; 
    mean_RH_NC_C_key, 2, 2];
err = [err_RH_NC_C_grp; err_RH_NC_C_key];

for i=1:size(xdata,1)
    h = bar(x1(i),xdata(i,1),0.5,'EdgeColor',[0 0 0],'LineWidth',1.5);
    if i == 1, hold on, end
    if xdata(i,3) == 1
        col = [82,107,218]/255;
    elseif xdata(i,3) == 2
        col = [248,117,99]/255;
    end
    set(h, 'FaceColor', col)
end
hold on;
ngroups = 1;
nbars = size(xdata, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    errorbar(x1(i), xdata(i), err(i), 'k', 'LineWidth',1.5,'linestyle', 'none');
end
legend(cndNames,'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('Noncrowding-crowding sensitivity difference')
xlabel('Condition')
ylabel('Sensitivity difference (d'')')
ylim([0 1.5])

%% PLOT D PRIMES / WITH GRASP / FOR LEFT HEMI
cndNames = {'NC, key','NC, grasp','C, key','C, grasp'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([LH_NC_dp_key; LH_NC_dp_grp; LH_C_dp_key; LH_C_dp_grp]); 
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

figure
x = [mean(LH_NC_dp_key); mean(LH_NC_dp_grp); mean(LH_C_dp_key); mean(LH_C_dp_grp)];
err = [std(LH_NC_dp_key)/sqrt(nSubj); std(LH_NC_dp_grp)/sqrt(nSubj); std(LH_C_dp_key)/sqrt(nSubj); std(LH_C_dp_grp)/sqrt(nSubj)]; 
h = bar(x); 
h(1).FaceColor = [135,205,215]/255;
hold on;
ngroups = 1;
nbars = size(x, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x1 = i;
    errorbar(x1, x(i), err(i), 'k', 'linestyle', 'none');
end
legend({'Average','SEM'},'AutoUpdate','off','location','northeast');

set(gca,gcaOpts{:})
title('Orientation sensitivity (d''), n=5')
xlabel('Condition')
ylabel('d-prime')
ylim([0 2.7])



%% PLOT D PRIME COMPARISON, GRASP VS NO GRASP / CATCH TRIALS
cndNames = {'NC','C'};
gcaOpts = {'XTick',1:length(cndNames),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',2,'fontsize',14};

figure
h = bar([LH_NC_diff; LH_C_diff]); 
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
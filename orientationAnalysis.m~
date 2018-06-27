%% DATASET INFO
% 1 Sydney
% 2 James
% 3 Christian
% 4 Bethany
% 5 Jianfei
% added baseline, increased maxori to 20, changed delta to 20%
% 6 Andrea
% 7 Ryan
% 8 Sydney
% 9 Bethany

% calculate w/in subjects left-right hemifield threshold difference 
% code T1 T2 - check jianfei's thing & Bekkering paper orientation change
% detection to avoid apparent motion 
% maybe reduce cntrast instead of mask? 

% move target such that its only in the upper hemifield at various radial
% positions

%%
clear all
close all
% parentDir = '~/Bethany/paclab';
parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));
dataDir = sprintf('%s/Subject_folders',parentDir);

nCnd = 8; % num staircases
subj = sort(strsplit(ls(dataDir)));
subj = subj(2:end-3); % ALL SUBJ
% subj = subj(2:5); % NO BASELINE
subj = subj(6:end); % WITH BASELINE

for s = 1:length(subj)
    curSubj = subj{s};
    file = mySubFiles(sprintf('%s/%s',dataDir,curSubj),'.m',18);
    data = load(file{:});
    
    stimRev = data.stimulusReversal;   
    nReverse(s,:) = data.nReverse;
    trialsPerStair = data.trial;
    acc = data.acc;
    acc = acc(1:max(trialsPerStair),:); % already split up by staircase
    
    for i = 1:nCnd % num staircases
        accVec = acc(1:trialsPerStair(i),i);
        accuracy(s,i) = sum(accVec)/trialsPerStair(i);
        
        finalRev(s,i) = stimRev(i,nReverse(s,i));
        final6Rev(:,i,s) = stimRev(i,nReverse(s,i)-5:nReverse(s,i));
        
        sumReversal(i) = sum(stimRev(i,nReverse(s,i)-5:nReverse(s,i)));
        stairmean(s,i) = sumReversal(i)/(nReverse(s,i)-3);
        StandardDev(s,i) = std(stimRev(i,4:nReverse(s,i)));
    end
    
    rspRatio = data.rspRatio;
    percRight(s) = rspRatio(2)/sum(rspRatio);
    
    % BIAS ANALYSIS SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trials = data.trials;
    rspKey = data.rspKey;
    rspKey = rspKey(1:trials);
    r1 = data.r1;
    r1 = r1(1:trials);
    hemiIdx = data.hemiIndex;
    hemiIdx = hemiIdx(1:trials);
    flankerIdx = data.flankerIndex;
    flankerIdx = flankerIdx(1:trials);
    
    rightTarget = find(r1 > 0);
    leftTarget = find(r1 < 0);
    rightPress = find(rspKey == 1);
    leftPress = find(rspKey == 0);
    rightHemi = find(hemiIdx == 1);
    leftHemi = find(hemiIdx == -1);
    rightFlanker = find(flankerIdx == 1);
    leftFlanker = find(flankerIdx == -1);
    
    % BIAS BASED ON HEMIFIELD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % for all trials in left hemi, how many - presses?
    counter = 0;
    for i = 1:length(leftHemi)
        if find(leftPress == leftHemi(i))
            counter = counter + 1;
        end
    end
    propLHemiLPress(s) = counter/length(leftHemi);
    
    % for all trials in left hemi, how many - targets?
    counter = 0;
    for i = 1:length(leftHemi)
        if find(leftTarget == leftHemi(i))
            counter = counter + 1;
        end
    end
    propLHemiLTarg(s) = counter/length(leftHemi);
    
    % for all trials in right hemi, how many + presses?
    counter = 0;
    for i = 1:length(rightHemi)
        if find(rightPress == rightHemi(i))
            counter = counter + 1;
        end
    end
    propRHemiRPress(s) = counter/length(rightHemi);
    
    % for all trials in right hemi, how many + targets?
    counter = 0;
    for i = 1:length(rightHemi)
        if find(rightTarget == rightHemi(i))
            counter = counter + 1;
        end
    end
    propRHemiRTarg(s) = counter/length(rightHemi);
    
    % BIAS BASED ON FLANKERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % for all trials with + flankers, how many responded +?
    counter = 0;
    for i = 1:length(rightFlanker)
        if find(rightPress == rightFlanker(i))
            counter = counter + 1;
        end
    end
    propRFlankerRPress(s) = counter/length(rightFlanker);
    
    % for all trials with + flankers, how many were +?
    counter = 0;
    for i = 1:length(rightFlanker)
        if find(rightTarget == rightFlanker(i))
            counter = counter + 1;
        end
    end
    propRFlankerRTarg(s) = counter/length(rightFlanker);
    
    % for all trials with - flankers, how many responded -?
    counter = 0;
    for i = 1:length(leftFlanker)
        if find(leftPress == leftFlanker(i))
            counter = counter + 1;
        end
    end
    propLFlankerLPress(s) = counter/length(leftFlanker);
    
    % for all trials with - flankers, how many were -?
    counter = 0;
    for i = 1:length(leftFlanker)
        if find(leftTarget == leftFlanker(i))
            counter = counter + 1;
        end
    end
    propLFlankerLTarg(s) = counter/length(leftFlanker);
    
    % for all trials with + flankers, how many responded -?
    counter = 0;
    for i = 1:length(rightFlanker)
        if find(leftPress == rightFlanker(i))
            counter = counter + 1;
        end
    end
    propRFlankerLPress(s) = counter/length(rightFlanker);
    
    % for all trials with + flankers, how many were -?
    counter = 0;
    for i = 1:length(rightFlanker)
        if find(leftTarget == rightFlanker(i))
            counter = counter + 1;
        end
    end
    propRFlankerLTarg(s) = counter/length(rightFlanker);
    
    % for all trials with - flankers, how many responded +?
    counter = 0;
    for i = 1:length(leftFlanker)
        if find(rightPress == leftFlanker(i))
            counter = counter + 1;
        end
    end
    propLFlankerRPress(s) = counter/length(leftFlanker);
    
    % for all trials with - flankers, how many were +?
    counter = 0;
    for i = 1:length(leftFlanker)
        if find(rightTarget == leftFlanker(i))
            counter = counter + 1;
        end
    end
    propLFlankerRTarg(s) = counter/length(leftFlanker);
end

grandMean = mean(stairmean,1); 
cndNames = {'Subj 1' 'Subj 2' 'Subj 3' 'Subj 4'};
gcaOpts = {'XTick',1:4,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

%% REORGANIZING & COLLAPSING DATA

updownMean = [mean(finalRev(:,1:2),2), mean(finalRev(:,3:4),2), mean(finalRev(:,5:6),2), mean(finalRev(:,7:8),2)];
sepCrwd = [updownMean(:,1) updownMean(:,3) updownMean(:,2) updownMean(:,4)];

%% PLOT: LEFT-RIGHT WITHIN SUBJECT COMPARISON

x = final6Rev(:,:,1);
y = LRdiff(:,:,1);

UDmean = [mean(final6Rev(:,1:2,:),2), mean(final6Rev(:,3:4,:),2), mean(final6Rev(:,5:6,:),2), mean(final6Rev(:,7:8,:),2)];
% 6 x 4 x 4: reversals x staircases x subjects

LRdiff2 = [UDmean(:,3,:)-UDmean(:,1,:), UDmean(:,4,:)-UDmean(:,2,:)];
% 6 x 2 x 4: reversals x staircases x subjects

meanLRdiff2 = mean(LRdiff2,1);

%%

% final6Rev: reversals x staircases x subjects
for i = 1:size(final6Rev,2)/2
    LRdiff(:,i,:) = final6Rev(:,i+4,:) - final6Rev(:,i,:);
end

withinSubjAvg = squeeze(mean(LRdiff,1))';
withinSubjSEM = squeeze(std(LRdiff,1)/sqrt(size(LRdiff,1)))';

for j = 1:size(LRdiff,3)
    UD(:,:,j) = [LRdiff(:,1,j), LRdiff(:,3,j); LRdiff(:,2,j), LRdiff(:,4,j)]; 
end
withinSubjUDAvg = squeeze(mean(UD,1))';
withinSubjUDSEM = squeeze(std(UD,1)/sqrt(size(UD,1)))';

%%
cndNames = {'D/NC','U/NC','D/C','U/C'};
gcaOpts = {'XTick',1:nCnd,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
bar(withinSubjAvg')
hold on;
ngroups = size(withinSubjAvg, 2);
nbars = size(withinSubjAvg, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - (groupwidth/2) + (2*i-1) * (groupwidth / (2*nbars));
    errorbar(x, withinSubjAvg(i,:), withinSubjSEM(i,:), 'k', 'linestyle', 'none');
end
legend(subj,'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Left-right hemifield bias comparisons, avged over last 6 reversals')
xlabel('Condition/Staircase')
ylabel('Threshold difference')
ylim([-10 10])

cndNames = {'NC','C'};
gcaOpts = {'XTick',1:nCnd,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
bar(withinSubjUDAvg')
hold on;
ngroups = size(withinSubjUDAvg, 2);
nbars = size(withinSubjUDAvg, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - (groupwidth/2) + (2*i-1) * (groupwidth / (2*nbars));
    errorbar(x, withinSubjUDAvg(i,:), withinSubjUDSEM(i,:), 'k', 'linestyle', 'none');
end
legend(subj,'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Left-right hemifield bias comparisons, avged over last 6 reversals')
xlabel('Condition/Staircase')
ylabel('Threshold difference')
ylim([-10 10])

%% PLOT THRESHOLDS
cndNames = {'L/NC','R/NC','L/C','R/C'};
gcaOpts = {'XTick',1:nCnd,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
bar(sepCrwd'); 
legend(subj,'AutoUpdate','off','location','northeast');
set(gca,gcaOpts{:})
title('Avg of last 6 reversal values, avged over stair dir')
xlabel('Condition/Staircase')
ylabel('Avg threshold (deg)')
ylim([0 20])

%% PLOT / for all trials w + flankers, how many - press?
figure
bar(propLTargetRFlanker)
set(gca,gcaOpts{:})
title('for all + flankers, % of - response')
xlabel('Subject')
ylabel('Percentage left (-) response')
ylim([0 1])
hold on
bar(propLPressRFlanker, 0.4, 'FaceAlpha',0.5)
legend({'correct' 'actual'},'AutoUpdate','off','location','northeast');
hold off

%% PLOT BIAS BASED ON HEMIFIELD
cndNames = {'Subj 1' 'Subj 2' 'Subj 3' 'Subj 4' 'Subj 5' 'Subj 6' 'Subj 7' 'Subj 8' 'Subj 9'};
gcaOpts = {'XTick',1:length(subj),'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
bar(propLHemiLTarg)
set(gca,gcaOpts{:})
title('% - (left) keypress with stimulus in left hemifield')
xlabel('Subject')
ylabel('Percentage - (left) response')
ylim([0 1])
hold on
bar(propLHemiLPress, 0.4, 'FaceAlpha',0.5)
legend({'prop of - targ in all L hemi' 'prop of - press in all L hemi'},'AutoUpdate','off','location','northeast');
hold off

figure
bar(propRHemiRTarg)
set(gca,gcaOpts{:})
title('% + (right) keypress with stimulus in right hemifield')
xlabel('Subject')
ylabel('Percentage - (left) response')
ylim([0 1])
hold on
bar(propRHemiRPress, 0.4, 'FaceAlpha',0.5)
legend({'prop of + targ in all R hemi' 'prop of + press in all R hemi'},'AutoUpdate','off','location','northeast');
hold off

%% BIAS BASED ON FLANKER

figure
bar(propRFlankerRTarg)
set(gca,gcaOpts{:})
title('% + (right) keypress in all trials with + flankers')
xlabel('Subject')
ylabel('Percentage + (right) response')
ylim([0 1])
hold on
bar(propRFlankerRPress, 0.4, 'FaceAlpha',0.5)
legend({'prop of + targ with + flanker' 'prop of + press with + flanker'},'AutoUpdate','off','location','northeast');
hold off

figure
bar(propLFlankerLTarg)
set(gca,gcaOpts{:})
title('% - (left) keypress in all trials with - flankers')
xlabel('Subject')
ylabel('Percentage - (left) response')
ylim([0 1])
hold on
bar(propLFlankerLPress, 0.4, 'FaceAlpha',0.5)
legend({'prop of - targ with - flanker' 'prop of - press with - flanker'},'AutoUpdate','off','location','northeast');
hold off

figure
bar(propRFlankerLTarg)
set(gca,gcaOpts{:})
title('% - (left) keypress in all trials with + flankers')
xlabel('Subject')
ylabel('Percentage - (left) response')
ylim([0 1])
hold on
bar(propRFlankerLPress, 0.4, 'FaceAlpha',0.5)
legend({'prop of - targ with + flanker' 'prop of - press with + flanker'},'AutoUpdate','off','location','northeast');
hold off

figure
bar(propLFlankerRTarg)
set(gca,gcaOpts{:})
title('% + (right) keypress in all trials with - flankers')
xlabel('Subject')
ylabel('Percentage + (right) response')
ylim([0 1])
hold on
bar(propLFlankerRPress, 0.4, 'FaceAlpha',0.5)
legend({'prop of + targ with - flanker' 'prop of + press with - flanker'},'AutoUpdate','off','location','northeast');
hold off

%% PLOTS / BAR OF AVG REVERSAL VALUE (THRESHOLD)
% cndNames = {'D/L/-','U/L/-','D/L/+','U/L/+','D/R/-','U/R/-','D/R/+','U/R/+'};
cndNames = {'D/L/NC','U/L/NC','D/L/C','U/L/C','D/R/NC','U/R/NC','D/R/C','U/R/C'};
gcaOpts = {'XTick',1:nCnd,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

figure
h = bar(stairmean');
h(1).FaceColor = [229,66,66]/255;
h(2).FaceColor = [35,169,181]/255;
h(3).FaceColor = [145 186 218]/255;
legend(subj,'AutoUpdate','off','location','northeast');
hold on;
ngroups = size(stairmean, 2);
nbars = size(stairmean, 1);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - (groupwidth/2) + (2*i-1) * (groupwidth / (2*nbars));
    errorbar(x, stairmean(i,:), StandardDev(i,:), 'k', 'linestyle', 'none');
end
set(gca,gcaOpts{:})
title('Avg of last 6 reversal values')
xlabel('Condition/Staircase')
ylabel('Avg threshold (deg)')
ylim([0 20])

cndNames = {'L/NC','R/NC','L/C','R/C'};
gcaOpts = {'XTick',1:nCnd,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};
figure
bar(sepCrwd'); 
set(gca,gcaOpts{:})
title('Avg of last 6 reversal values, avged over subjects')
xlabel('Condition/Staircase')
ylabel('Avg threshold (deg)')
ylim([0 20])

%% PLOTS / BAR OF FINAL REVERSAL VALUE (THRESHOLD)
% cndNames = {'D/L/-','U/L/-','D/L/+','U/L/+','D/R/-','U/R/-','D/R/+','U/R/+'};
% cndNames = {'D/L/NC','U/L/NC','D/L/C','U/L/C','D/R/NC','U/R/NC','D/R/C','U/R/C'};
% gcaOpts = {'XTick',1:nCnd,'XTickLabel',cndNames,'box',...
%     'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};
% 
% figure
% h = bar(finalRev');
% h(1).FaceColor = [229,66,66]/255;
% h(2).FaceColor = [35,169,181]/255;
% h(3).FaceColor = [145 186 218]/255;
% legend(subj,'AutoUpdate','off','location','northeast');
% set(gca,gcaOpts{:})
% title('Final reversal values, excluding first 3')
% xlabel('Condition/Staircase')
% ylabel('Avg threshold (deg)')
% ylim([0 20])

%% TTEST TO TEST + vs - FLANKERS
clear h
idx = 0;
for i = [1,2,5,6]
    idx = idx+1;
    cnd1 = stairmean(:,i);
    cnd2 = stairmean(:,i+2);
    [h(idx),p(idx)] = ttest(cnd1,cnd2);
end

%% COMPARING STAIRCASES
subjStairMean = mean(stairmean, 1);
subjStairDev = std(stairmean,1); 

bar(subjStairMean)
hold on
errorbar(1:8,subjStairMean,subjStairDev,'.')

%%
figure
hold on
for i = 2:size(finalRev,1)
    plot(finalRev(i,:),'o')
end

% available variables: 
% trials
% nReverse
% trial
% stimulusReversal
% rspKey
% hemiIndex
% flankerIndex
% acc
% stairOrder

% TO DO: cross-correlate stairOrder and acc, hemiIndex and acc,
% flankerIndex and acc, hemiIndex and flankerIndex and acc to spot trends

% PUT BACK IN WHEN BLOCKS ARE A THING
%     blocks = mySubFiles(sprintf('%s/%s',dataDir,curSubj),'.m',18);
%     for i = 1:length(blocks)
%         data = load(blocks{i});
%         stimRev = data.stimulusReversal;
%         for j = 1:size(stimRev,1)
%             finalRev(j) = stimRev(j,max(find(stimRev(j,:) ~= 0)));
%         end
%         avgRev(i) = mean(finalRev);
%     end
%     allAvgRev(s,:) = avgRev;
%% DEVELOPMENT

for i = 1:8 % num staircases
    accVec = acc(1:trialsPerStair(i),i);
    accuracy(i) = sum(accVec)/trialsPerStair(i);
end
%%
i = 1;
data = load(blocks{i});
stimRev = data.stimulusReversal;
for j = 1:size(stimRev,1)
    finalRev(j) = stimRev(j,max(find(stimRev(j,:) ~= 0)));
end
avgRev = mean(finalRev);

rspKey = data.rspKey;
trials = data.trials;
rspKey = rspKey(1:trials);


%%
addpath(genpath(parentDir));
results = load(sprintf('%s/%d/%s',dataDir,1,blocks{1}));
% block2 = load(sprintf('%s/Subject_folders/1_block2/threshold.mat',parentDir));

stimulusReversal = results.stimulusReversal;
nReverse = results.nReverse;

plot(stimulusReversal(1,1:nReverse(1)));hold on;
plot(stimulusReversal(2,1:nReverse(2)));hold on;
%%

for i = 1:8
    sumReversal(i) = sum(stimRev(i,4:nReverse(i)));
    stairmean(i) = sumReversal(i)/(nReverse(i)-3);
    StandardDev(i) = std(stimRev(i,4:nReverse(i)));
end

plot(stimRev(1,1:nReverse(1)));hold on;
plot(stimRev(2,1:nReverse(2)));hold on;

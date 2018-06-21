%% DATASET INFO
% 1 = Sydney
% 2 = James
% 3 = Christian

%%
clear all
close all
% parentDir = '~/Bethany/paclab';
parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));

nCnd = 8; % num staircases

dataDir = sprintf('%s/Subject_folders',parentDir);
subj = sort(strsplit(ls(dataDir)));
% subj = subj(2:end-3);
subj = subj(3:end-3); % excluding sydney
for s = 1:length(subj)
    curSubj = subj{s};
    file = mySubFiles(sprintf('%s/%s',dataDir,curSubj),'.m',18);
    data = load(file{:});
    stimRev = data.stimulusReversal;
    for j = 1:size(stimRev,1)
        finalRev(s,j) = stimRev(j,max(find(stimRev(j,:) ~= 0)));
    end
    
    nReverse(s,:) = data.nReverse;

    trialsPerStair = data.trial;
    acc = data.acc;
    acc = acc(1:max(trialsPerStair),:); % already split up by staircase
    
    for i = 1:nCnd % num staircases
        accVec = acc(1:trialsPerStair(i),i);
        accuracy(s,i) = sum(accVec)/trialsPerStair(i);
    end
    for i = 1:nCnd
        sumReversal(i) = sum(stimRev(i,4:nReverse(s,i)));
        stairmean(s,i) = sumReversal(i)/(nReverse(s,i)-3);
        StandardDev(s,i) = std(stimRev(i,4:nReverse(s,i)));
%         allRevData(s,i) = stimRev(i,4:min(nReverse));
    end
    
    rspRatio = data.rspRatio;
    percRight(s) = rspRatio(2)/sum(rspRatio);
    
    % analyzing hemifield dependence of response bias
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
    
    % for all + presses, how many occurred in the right hemi? 
    counter = 0;
    for i = 1:length(rightPress)
        if find(rightHemi == rightPress(i))
            counter = counter + 1;
        end
    end
    propRPressRHemi(s) = counter/length(rightPress);
    
    % for all + targets, how many occurred in the right hemi? 
    counter = 0;
    for i = 1:length(rightTarget)
        if find(rightHemi == rightTarget(i))
            counter = counter + 1;
        end
    end
    propRTargRHemi(s) = counter/length(rightTarget);
    
%     % for all - presses, how many occurred in the left hemi? 
%     counter = 0;
%     for i = 1:length(leftPress)
%         if find(leftHemi == leftPress(i))
%             counter = counter + 1;
%         end
%     end
%     propLPressLHemi(s) = counter/length(leftPress);
%     
%     % for all - targets, how many occurred in the left hemi? 
%     counter = 0;
%     for i = 1:length(leftTarget)
%         if find(leftHemi == leftTarget(i))
%             counter = counter + 1;
%         end
%     end
%     propLTargLHemi(s) = counter/length(leftTarget);
    
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
    
    % for all trials with + flankers, how many were +?
    counter = 0;
    for i = 1:length(leftFlanker)
        if find(leftTarget == leftFlanker(i))
            counter = counter + 1;
        end
    end
    propLFlankerLTarg(s) = counter/length(leftFlanker);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     % for all + presses, how many had + flankers?? / CORRECT
%     counter = 0;
%     for i = 1:length(rightPress)
%         if find(rightFlanker == rightPress(i))
%             counter = counter + 1;
%         end
%     end
%     propRPressRFlanker(s) = counter/length(rightPress);
%     
%     % for all + targets, how many had + flankers? 
%     counter = 0;
%     for i = 1:length(rightTarget)
%         if find(rightFlanker == rightTarget(i))
%             counter = counter + 1;
%         end
%     end
%     propRTargAndRFlank(s) = counter/length(rightTarget);
%     
%     % for all trials with - presses, how many + flankers? 
% %     counter = 0;
% %     for i = 1:length(leftPress)
% %         if find(rightFlanker == leftPress(i))
% %             counter = counter + 1;
% %         end
% %     end
% %     propRFlankerLPress(s) = counter/length(leftPress);
%     
%     % for all trials w/ + flankers, how many - presses? 
%     counter = 0;
%     for i = 1:length(rightFlanker)
%         if find(rightFlanker(i) == leftPress)
%             counter = counter + 1;
%         end
%     end
%     propLPressRFlanker(s) = counter/length(rightFlanker);
%     
%     % now check actual frequency of - targets during all trials w/ + flankers
%     counter = 0;
%     for i = 1:length(rightFlanker)
%         if find(rightFlanker(i) == leftTarget)
%             counter = counter + 1;
%         end
%     end
%     propLTargetRFlanker(s) = counter/length(rightFlanker);
%     
%     % how many - targets occurred with + flankers? to check accuracy of the + presses during + flankers
%     counter = 0;
%     for i = 1:length(leftTarget)
%         if find(rightFlanker == leftTarget(i))
%             counter = counter + 1;
%         end
%     end
%     propLTargAndRFlank(s) = counter/length(leftTarget);
end

cndNames = {'Subj 1' 'Subj 2' 'Subj 3' 'Subj 4'};
gcaOpts = {'XTick',1:4,'XTickLabel',cndNames,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

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

%% PLOT / BAR OF LOCAL ORI & GLOBAL HEMI BIAS

% cndNames = {'Subj 1' 'Subj 2' 'Subj 3' 'Subj 4'};
gcaOpts = {'XTick',1:4,'XTickLabel',subj,'box',...
    'off','tickdir','out','fontname','Helvetica','linewidth',1.5,'fontsize',14};

% figure
% bar(propRTargRHemi)
% set(gca,gcaOpts{:})
% title('% + (right) keypress with stimulus in right hemifield')
% xlabel('Subject')
% ylabel('Percentage + (right) response')
% ylim([0 1])
% hold on
% bar(propRPressRHemi, 0.4, 'FaceAlpha',0.5)
% legend({'prop of + targ in R hemi' 'prop of + press in R hemi'},'AutoUpdate','off','location','northeast');
% hold off
% 
% 
% figure
% bar(propLTargLHemi)
% set(gca,gcaOpts{:})
% title('% - (left) keypress with stimulus in left hemifield')
% xlabel('Subject')
% ylabel('Percentage - (left) response')
% ylim([0 1])
% hold on
% bar(propLPressLHemi, 0.4, 'FaceAlpha',0.5)
% legend({'prop of - targ in L hemi' 'prop of - press in L hemi'},'AutoUpdate','off','location','northeast');
% hold off

% BIAS BASED ON HEMIFIELD
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

%%
x = 0.5*1.3
x = x*1.3

x = 20-20*.3
x = x-x*.3

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
bar(propRFlankerRPress, 0.4, 'FaceAlpha',0.5)
legend({'prop of - targ with - flanker' 'prop of - press with - flanker'},'AutoUpdate','off','location','northeast');
hold off

%% PLOTS / BAR OF AVG REVERSAL VALUE (THRESHOLD)
cndNames = {'D/L/-','U/L/-','D/L/+','U/L/+','D/R/-','U/R/-','D/R/+','U/R/+'};
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
title('Avg reversal values, excluding first 3')
xlabel('Condition/Staircase')
ylabel('Avg threshold (deg)')
ylim([0 12])

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

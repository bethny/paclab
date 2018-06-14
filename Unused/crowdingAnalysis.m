%% LOAD & SORT ANSWERS
clear all
close all

atLab = 1;
parentDir = {'~/code/pac','~/Bethany/paclab'};
dataDir = {'~/code/pac/Data/Psychophysics', '~/Bethany/paclab/Data/Psychophysics'};
addpath(genpath(parentDir{atLab+1}));

%%
subj = ls(dataDir{atLab+1});
subj = strsplit(subj);
for s = 1:length(subj)-1
    curSubj = subj{s};
    blocks = mySubFiles(sprintf('%s/%s',dataDir{atLab+1},curSubj),curSubj,1);
    trialInfo(:,:,s) = getTrialInfo(blocks); % crowding str, angle, accuracy % FIX bc not all subj have equal numbers of blocks
end

% sort out the various crowding strengths
noCrwd = trialInfo(find(trialInfo(:,1,:) == 0),:,:);
lowCrwd = trialInfo(find(trialInfo(:,1,:) == 1),:,:);
strongCrwd = trialInfo(find(trialInfo(:,1,:) == 2),:,:);

% calculate perc of correct rsp per cond
rspNC = calcCor(noCrwd); % no crowding
rspLC = calcCor(lowCrwd); % low crowding
rspSC = calcCor(strongCrwd); % strong crowding
% [output] = gaborDisplay;

%% LOAD & SORT ANSWERS
clear all
close all

codeDir = '~/code';
parentDir = '~/code/pac/Data/Psychophysics';
addpath(genpath(parentDir))
addpath(genpath(codeDir))
%%
subj = ls(parentDir);
subj = strsplit(subj);
for s = 1:length(subj)-1
    curSubj = subj{s};
    blocks = mySubFiles(sprintf('%s/%s',parentDir,curSubj),curSubj,1);
    trialInfo(:,:,s) = getTrialInfo(blocks); % crowding str, angle, accuracy
end

% sort out the various crowding strengths
noCrwd = trialInfo(find(trialInfo(:,1,:) == 0),:,:);
lowCrwd = trialInfo(find(trialInfo(:,1,:) == 1),:,:);
strongCrwd = trialInfo(find(trialInfo(:,1,:) == 2),:,:);

% calculate perc of correct rsp per cond
rspNC = calcCor(noCrwd); % no crowding
rspLC = calcCor(lowCrwd); % low crowding
rspSC = calcCor(strongCrwd); % strong crowding
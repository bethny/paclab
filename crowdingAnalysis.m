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
    trialInfo(:,:,s) = getTrialInfo(blocks);
end

noCrowding = trialInfo(find(trialInfo(:,1,:) == 0),:,:);
lowCrowding = trialInfo(find(trialInfo(:,1,:) == 1),:,:);
strongCrowding = trialInfo(find(trialInfo(:,1,:) == 2),:,:);

rspNC = calcCor(noCrowding);
rspLC = calcCor(lowCrowding);
rspSC = calcCor(strongCrowding);
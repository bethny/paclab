function trialInfo = getTrialInfo(blocks)

% INPUT: .mat file containing block data (response, RT, stimset)
% OUTPUT: 

nBlocks = length(blocks);
subjRT = [];
subjRsp = [];
subjStim = [];
for b = 1:nBlocks
    curFile = load(blocks{b});
    subjRT = [subjRT; curFile.output.rsp.RT'];
    subjRsp = [subjRsp; curFile.output.rsp.keyName'];
    curStim = curFile.output.cnd(1:length(curFile.output.rsp.RT'),:); 
    subjStim = [subjStim; curStim];
end

% convert key codes to numbers
numRsp = zeros(size(subjRsp));
numRsp(strcmp(subjRsp,'LeftArrow')) = 1;
numRsp(strcmp(subjRsp,'RightArrow')) = 2;
numRsp(strcmp(subjRsp,'none')) = 0;

% combine into trialInfo array: cols = crowding str, angle, resp, correctness
trialInfo = cat(2,subjStim,numRsp,zeros(size(numRsp)));
for i = 1:length(trialInfo)
    if trialInfo(i,2) < 0 && trialInfo(i,4) == 1 || trialInfo(i,2) > 0 && trialInfo(i,4) == 2
        trialInfo(i,5) = 1;
    else
        trialInfo(i,5) = 0;
    end
end

% final output: cols = crowding str, angle, hemifield, correctnes
trialInfo = trialInfo(:,[1:3 5]);

end
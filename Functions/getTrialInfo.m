function trialInfo = getTrialInfo(blocks)

nBlocks = length(blocks);
for b = 1:nBlocks
    curFile = load(blocks{b});
    subjRt(b,:) = curFile.output.rsp.RT;
    subjRsp(b,:) = curFile.output.rsp.keyName;
    subjStim(:,:,b) = curFile.output.cnd; % 144 x 2 x 2        
end
relStim = subjStim(1:size(subjRt,2),:,:); % positive = left, neg = right  
ugh = permute(relStim,[1 3 2]);
vecRelStim = reshape(ugh,[],size(relStim,2),1);

vecRsp = reshape(subjRsp,1,[])';
numRsp = zeros(size(vecRsp));
numRsp(find(strcmp(vecRsp,'LeftArrow'))) = 1;
numRsp(find(strcmp(vecRsp,'RightArrow'))) = 2;
numRsp(find(strcmp(vecRsp,'none'))) = 0;

trialInfo = cat(2,vecRelStim,numRsp,zeros(size(numRsp)));
leftIdx = find(trialInfo(:,2) > 0);
rightIdx = find(trialInfo(:,2) < 0);
for l = 1:length(leftIdx)
    if trialInfo(leftIdx(l),3) == 1
        trialInfo(leftIdx(l),4) = 1;
    end
end
for r = 1:length(rightIdx)
    if trialInfo(rightIdx(r),3) == 2
        trialInfo(rightIdx(r),4) = 1;
    end
end  

trialInfo = trialInfo(:,[1 2 4]);

end
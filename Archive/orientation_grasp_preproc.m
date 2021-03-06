% tracker records at 160 Hz

% 1/23/12 UPDATED 7/2016
% By Jeff Moher
% Program for analyzing reaching movement in a singleton detection paradigm
% Set custom parameters at beginning of code, and also search for "SPECIAL
% SECTION" headings for sections that have to be updated depending on
% specifics of a given project
% Updated by Jianfei Guo 11/25/16
clear all;
warning off;

parentDir = '~/code/pac/paclab';
addpath(genpath(parentDir));
dataDir = sprintf('%s/Subject_folders',parentDir);
subj = sort(strsplit(ls(dataDir)));
subj = subj(end);
nSubj = length(subj);

% Instructions.  If examining images and adjusting thresholds by hand, use
% the following keys:
% 's' - stay at this image
% 'f' - foward to following images
% 'z' - skip to next image, current is OK
% 'b' - go back one image
% 'd' - delete the current image, leave a note as to why.  will be saved to file

% Note that this assumes a particular file structure - in a folder called
% "Data", there should be one folder for each participant.  In each of
% these folders, there should be two files for each block - one with
% reaching data, one with other variables.

% Note - keyboard output to the command line is suppressed, so if it
% crashes, you may need to press ctrl+c to get the keyboard back.

% Enter custom variables here:

% This is where you specify which variables you want to extract from your
% non-reaching data file.  For example, let's say you need two variables -
% trial number and distractor presence.   You would enter them here (order
% doesn't matter) as you see below, where the '.name' specifies what you
% want to call that variable, and the '.number' specifies which column has
% that variable's data.

allVars = struct;
allVars(1).name = 'targetID';
allVars(1).number = '3';
allVars(2).name = 'baseOri';
allVars(2).number = '4';
allVars(3).name = 'crowding';
allVars(3).number = '5';
allVars(4).name = 'reachedTo';
allVars(4).number = '6';
allVars(5).name = 'timeElapsed';
allVars(5).number = '9';
allVars(6).name = 'angle';
allVars(6).number = '11';
allVars(7).name = 'acc';
allVars(7).number = '12';
allVars(8).name = 'acc_grasp';
allVars(8).number = '13';
allVars(9).name = 'markerWait';
allVars(9).number = '14';

for s = 1:nSubj
    fprintf(sprintf('Now processing subject %d!\n',s));
    curSubj = subj{s};
    
    subjDir = sprintf('%s/%s',dataDir,curSubj);
    cd(subjDir);
    subjNum = str2num(curSubj(2:3));
    
    % Note that this just assigns numbers to the subject folders in
    % alphabetical order, so if you're first subject is #305, that will be
    % subject #1 for these purposes.
    % Various flags
    drawImages = 1; %1 if you want to see each threshold
    saveFiles = 1; %1 if you want to save the new files
    reUseThresholds = 0; %1 if you want to re-use thresholds from a previous pre-processing
    startMiddle = 0;% 1 if you want to start in the middle - this is uesful if you crash
    %in the middle of pre-processing, so you don't have to start from the
    %beginning.  It will prompt you asking which block and trial number you
    %want to start at
    drawInPixels = 1;% currently only works as 1
    originalTimeLimit = 1; % How long should each image be displayed before auto-advancing
    timeLimit = originalTimeLimit;
    % Specify the length of each block in the current study
    blockLength = 64;
    if subjNum == 5
        nBlocks = 2;
    else
        nBlocks = 4;
    end
    totalTrials = blockLength * nBlocks;
    % mark as "1" if you want to upsample data
    upSample = 0;
    samplingRate = 160; % orig 200
    
    % Current subject for analysis, if running only one subject at a time.
    %     for currentSub = curSubj
    %currentSub = 14;
    
    if startMiddle
        blockEnd = str2num(input('\n Enter last Block Number: ', 's'));
        trialEnd = str2num(input('\n Enter last Trial Number: ', 's'));
    end
    
    % Name of folder where you'll keep change files
    if ~isdir('ManualChanges')
        mkdir('ManualChanges');
    end
    subFolder = 'ManualChanges';
    
    ListenChar(2); %suppress keyboard output in the command line
    
    if saveFiles && ~reUseThresholds
        % % Create name for output file for Redo Trials (trials where visual
        % inspection suggests something is wrong) and open it
        dataFilename = [subFolder,filesep,'RedoTrials_',curSubj,'.txt'];
        dataFilepointer = fopen(dataFilename, 'w');
        %fprintf(dataFilepointer, '%s \t  %s \t %s \t  %s \n','Subject','Block', 'Trial', 'Notes');
        
        % Create name for output file for modified threshold Trials (trials where the movement thresholds were manually
        % changed) and open it
        dataFilename2 = [subFolder,filesep,'ThresholdAdjustments_',curSubj,'.txt'];
        dataFilepointer2 = fopen(dataFilename2, 'w');
        %fprintf(dataFilepointer2, '%s \t  %s \t %s \t %s \t %s \n','Subject','Block', 'Trial', 'NewStartIndex', 'NewEndIndex');
    elseif reUseThresholds %&& startMiddle
        % % Create name for output file for Redo Trials (trials where visual
        % inspection suggests something is wrong) and open it
        dataFilename = [subFolder,filesep,'RedoTrials_',curSubj,'.txt'];
        dataFilepointer = fopen(dataFilename, 'a');
        %fprintf(dataFilepointer, '%s \t  %s \t %s \t  %s \n','Subject','Block', 'Trial', 'Notes');
        
        % Create name for output file for modified threshold Trials (trials where the movement thresholds were manually changed)
        % and open it
        dataFilename2 = [subFolder,filesep,'ThresholdAdjustments_',curSubj,'.txt'];
        dataFilepointer2 = fopen(dataFilename2, 'a');
        %fprintf(dataFilepointer2, '%s \t  %s \t %s \t %s \t %s \n','Subject','Block', 'Trial', 'NewStartIndex', 'NewEndIndex');
    end
    
    %Change to path of data
    path = pwd;
    tempFiles = dir(path);
    dirFlags = [tempFiles.isdir];
    folders = tempFiles(dirFlags);
    folders(1:2) = [];
    
    if drawImages
        % Initiate two figures for examination - top left = movement itself, bottom
        % left = velocity profile
        scrsz = get(0,'ScreenSize');
        h = figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
        hh = figure('Position',[1 scrsz(2) scrsz(3)/2 scrsz(4)/4]);
        grid ON;
    end
    
    % For loop to analyze data for each subject
    for allFiles = subjNum % Can be changed to pre-process multiple subjects at once. Should use only if re-doing analysis
        % with pre-processed files
        % Code to call up files with info about which manual changes have been made
        if reUseThresholds % && startMiddle
            changeFile = [subFolder,filesep,'ThresholdAdjustments_',curSubj,'.txt'];
            discardFile = [subFolder,filesep,'RedoTrials_',curSubj,'.txt'];
            dataMatrix = load(changeFile);
            dataMatrix2 = load(discardFile);
            if ~isempty(dataMatrix)
                changesBlock = dataMatrix(:,2);
                changesTrial = dataMatrix(:,3);
                newStartIndices = dataMatrix(:,4);
                newEndIndices = dataMatrix(:,5);
            end
            
            if ~isempty(dataMatrix2)
                discardBlock = dataMatrix2(:,2);
                discardTrial = dataMatrix2(:,3);
                dropsCounter = 1;
            end
            
            changesCounter = 1;
            showTrial = 0;
        end
        
        clickedThreshold = 0;
        
        % within the loop, open each data folder in order to analyze each individual file
        cd(subjDir);
        files = dir(subjDir);
        files(1:2) = [];
        fileIndexStart = 1;
        if strcmp(files(1).name,'.DS_Store')
            files(1) = [];
        end
        
        % Which column in the data file is used to measure x & y coordinates
        % Currently in space (and thus inches - need to convert to cm)
        trialNumberMarker = 1;
        blockNumberMarker = 2;
        timeMarker = 2+1;
        x1Marker = 7+1;
        z1Marker = 9+1;
        y1Marker = 8+1;
        x2Marker = 10+1;
        z2Marker = 12+1;
        y2Marker = 11+1;
        sampleMarker = 13+1;
        
        % This is where the threshold is set for what determines the start and end of a movement. The min samples
        % portions determine how many consectuive samples must exceed (or fall below) the threshold to determine the
        % start/end of a movement.
        standardThreshold = 15; % cm/s
        standardLowerThreshold = 15;
        standardMinSamples = 1;
        standardMinEndSamples = 1;
        
        % Set minimum thresholds (in cm/s) to define the start and end of arm movement, i.e. what velocity (or acceleration,
        % if needed) must be exceeded to consider something an arm movement
        threshold = standardThreshold;
        lowerThreshold = standardLowerThreshold;
        
        % Set minimun number of consecutive samples over/under threshold required to accept the start/end of a movement
        minSamples = standardMinSamples;
        minEndSamples = standardMinEndSamples;
        
        % preallocate variables that will be recorded as totals for all blocks
        maxDeviation = zeros(totalTrials,1);
        maxDeviationUncorrected = zeros(totalTrials,1);
        maxAngularDeviation = zeros(totalTrials,1);
        deviationAtPeakVelocity = zeros(totalTrials,1);
        deviationAtPeakAcceleration = zeros(totalTrials,1);
        deviationAtPeakDeceleration = zeros(totalTrials,1);
        angleAtPeakVelocity = zeros(totalTrials,1);
        peakVelocity = zeros(totalTrials,1);
        peakAcceleration = zeros(totalTrials,1);
        turnAroundPoint = zeros(totalTrials,1);
        turnAroundRT = zeros(totalTrials,1);
        signOfCurvature = zeros(totalTrials,1);
        MT = zeros(totalTrials,1);
        RT = zeros(totalTrials,1);
        acc = zeros(totalTrials,1);
        noMovement = zeros(totalTrials,1);
        droppedTrials = zeros(totalTrials,1);
        angularXDeviation = zeros(totalTrials,1);
        rightwardX = zeros(totalTrials,1);
        signOfCurvatureX = zeros(totalTrials,1);
        overallSRate = [];
        AllXMovementPoints = zeros(totalTrials,1000);
        AllYMovementPoints = zeros(totalTrials,1000);
        
        % Preallocate variables that will be recorded separately for each block
        reachStartIndex = zeros(blockLength,1);
        reachEndIndex = zeros(blockLength,1);
        maxDevPoint = zeros(blockLength,1);
        maxCurrentX = zeros(blockLength,1);
        maxCurrentY = zeros(blockLength,1);
        startx = zeros(blockLength,1);
        starty = zeros(blockLength,1);
        endx = zeros(blockLength,1);
        endy = zeros(blockLength,1);
        peakDeceleration = zeros(blockLength,1);
        peakVelocityIndex = zeros(blockLength,1);
        peakAccelerationIndex = zeros(blockLength,1);
        peakDecelerationIndex = zeros(blockLength,1);
        
        % for loop to analyze individual data files for current subject
        for b = 1:nBlocks
            allVelocity = [];
            
            % Open two files for each run through the loop  - simple data file = file, and movement data file = file2
            data = [];
            movementData = [];
            
            file = [subjDir,filesep,num2str(subjNum),'_',num2str(b),'_graspExp_info.txt'];
            data2 = load(file);
            data = [data; data2];
            
            file2 = [subjDir,filesep,num2str(subjNum),'_',num2str(b),'_graspExp.txt']; % why doesn't this loop over blocks?
            data3 = load(file2);
            movementData = [movementData; data3];
            
            % delete crashed trials
            index = find(diff(data(:,1)) == 0);
            trial = data(index,1);
            block = data(index,2);
            data(index,:) = [];
            
            if upSample
                tempa = [1;diff(movementData(:,2))];
                temptrialNumbers = movementData(:,2);
                temptrialStartIndex = find(tempa);
                temptrialEndIndex = [temptrialStartIndex(2:blockLength) - 1;length(movementData)];
                
                for j = 1:length(temptrialStartIndex)
                    timeElapsed = movementData(temptrialEndIndex(j),timeMarker) - movementData(temptrialStartIndex(j),...
                        timeMarker);
                    upsampleRate = round(timeElapsed * samplingRate);
                    numSamples = temptrialEndIndex(j) - temptrialStartIndex(j);
                    if samplingRate > (numSamples/timeElapsed)
                        tempData = zeros(upsampleRate,9);
                        tempData(:,1) = movementData(temptrialStartIndex(j),1);
                        tempData(:,2) = movementData(temptrialStartIndex(j),2);
                        newpath=pwd;
                        cd(oripath);
                        tempData(:,3) = EZResample(movementData(temptrialStartIndex(j):temptrialEndIndex(j),3),upsampleRate,0);
                        tempData(:,4) = EZResample(movementData(temptrialStartIndex(j):temptrialEndIndex(j),4),upsampleRate,0);
                        tempData(:,5) = EZResample(movementData(temptrialStartIndex(j):temptrialEndIndex(j),5),upsampleRate,0);
                        tempData(:,6) = EZResample(movementData(temptrialStartIndex(j):temptrialEndIndex(j),6),upsampleRate,0);
                        tempData(:,7) = EZResample(movementData(temptrialStartIndex(j):temptrialEndIndex(j),7),upsampleRate,0);
                        tempData(:,8) = EZResample(movementData(temptrialStartIndex(j):temptrialEndIndex(j),8),upsampleRate,0);
                        tempData(:,9) = EZResample(movementData(temptrialStartIndex(j):temptrialEndIndex(j),9),upsampleRate,0);
                        movementData = [movementData;tempData];
                        cd(newpath);
                        %append to movement data
                    else
                        'help!!'
                        j
                        return;
                        %append current to movement data
                    end
                end
            else
                movementData = load(file2);
            end
            
            % Convert movement data to cm
            movementData(:,x1Marker) = movementData(:,x1Marker) * 2.54;
            movementData(:,y1Marker) = movementData(:,y1Marker) * 2.54;
            movementData(:,z1Marker) = movementData(:,z1Marker) * 2.54;
            movementData(:,x2Marker) = movementData(:,x2Marker) * 2.54;
            movementData(:,y2Marker) = movementData(:,y2Marker) * 2.54;
            movementData(:,z2Marker) = movementData(:,z2Marker) * 2.54;
            
            % delete crashed trials
            delete_index = [];
            ii = 1;
            for jj = 1:length(block)
                time = 0;
                flag = 0;
                while ~flag
                    if any(movementData(ii,2)==block(jj)) && any(movementData(ii,1)==trial(jj)) && (movementData(ii,3)>time)
                        delete_index=[delete_index,ii];
                        time = movementData(ii,3);
                    elseif any(movementData(ii,2)==block(jj))&& any(movementData(ii,1)==trial(jj))&& (movementData(ii,3)<time)
                        flag=1;
                    end
                    ii = ii+1;
                end
            end
            movementData(delete_index,:)=[];
            
            for j = 1:length(allVars)
                eval([allVars(j).name,' = data(:,',allVars(j).number,');']);
            end
            
            % Calculate speed for every sample in the file. A zero is added at the beginning because point 1 cm the file
            % will not have a speed, and the diff command will output the first speed
            % value, from point 2, as the first value in the scalar.  Adding a
            % zero will line up the vector correctly
            timeChange = diff(movementData(:,timeMarker));
            xArmChange = -diff(movementData(:,x1Marker)); % Distance between each sample in X dimension
            yArmChange = -diff(movementData(:,y1Marker)); % Distance between each sample in Y dimension
            zArmChange = -diff(movementData(:,z1Marker)); % Distance between each sample in Z dimension
            totalArmChange = sqrt(xArmChange.^2 + yArmChange.^2 + zArmChange.^2); % Calculate the 2D distance between each movement
            tempAllVelocity = [0;totalArmChange./timeChange]; % Calculate velocity by dividing distance by time change
            
            % Calculate acceleration for every sample in the file.  Add a zero at the
            % beginning for same logic above - first data point will not have
            % an acceleration
            velocityChange = diff(tempAllVelocity); % Change in velocity between each sample
            allAcceleration = [0;velocityChange./timeChange]; % Calculate acceleration by dividing velocity change by time change
            
            a = [1;diff(movementData(:,trialNumberMarker))];
            trialNumbers = movementData(:,trialNumberMarker);
            
            % Look for trial numbers - NOTE: this assings a trial start and
            % end index of "1" if it fails to find an instance of a given
            % trial number
            trialStartIndex = find(a);
            trialEndIndex = [trialStartIndex(2:end)  - 1;length(movementData)];
            
            % Reset velocities and accelerations calculated on the first sample
            % of each trial to zero, because they don't reflect a real change
            % in velocity.
            tempAllVelocity(trialStartIndex) = 0;
            allAcceleration(trialStartIndex) = 0;
            
            
            % Filter the velocity data for smoothness using a butterworth
            % filter - parameters specified within
            for xx = 1:length(trialStartIndex)
                % Parameters for Butterworth filter
                % Sampling rate for the reach tracker
                timeElapsedForThisTrial = movementData(trialEndIndex(xx),timeMarker) - movementData(trialStartIndex(xx),...
                    timeMarker);
                samplesForThisTrial = trialEndIndex(xx) - trialStartIndex(xx);
                sRate = round(samplesForThisTrial/timeElapsedForThisTrial);
                overallSRate = [overallSRate;sRate];
                % Half of sampling rate, for input to the filter
                halfSRate = sRate/2;
                % Frequency in HZ for high cut for Butterworth filter
                freqCut = 10;
                % "Window" that needs to be specified for Butterworth filter
                butterWindow = freqCut/halfSRate;
                % Nth order for Butterworth filter
                nthOrder = 2;
                if sRate>0 && samplesForThisTrial >= 6
                    [bttr,a] = butter(nthOrder,butterWindow);
                    tempx = filtfilt(bttr,a,tempAllVelocity(trialStartIndex(xx):trialEndIndex(xx)));
                    allVelocity(trialStartIndex(xx):trialEndIndex(xx)) = tempx;
                end
            end
            
            % Calculate timepoint of beginning of each trial - should be very
            % near zero
            trialStartTime = movementData(trialStartIndex,timeMarker);
            
            % Calculate where throughout the file velocity exceeds threshold
            velThresholdExceeded = allVelocity>threshold;
            
            % Calculate where velocity is exceeded for minimum # samples
            minSamplesCheck = ones(1,minSamples);
            minSamplesExceeded = strfind((velThresholdExceeded), minSamplesCheck);
            
            % Calculate where in the file velocity goes below threshold
            velLowerThresholdBelow = allVelocity<lowerThreshold;
            
            % Calculate where velocity goes below min threshold for minimum #
            % samples
            minEndSamplesCheck = ones(1,minEndSamples);
            minEndSamplesExceeded = strfind((velLowerThresholdBelow), minEndSamplesCheck);
            
            % In FOR loop for each trial, find the first place when velocity threshold was exceeded
            %  (if it was), and when it went back below threshold (if it did).
            %  Then, calculate curvature for each trial
            j = 1;
            
            while j <= (length(trialStartIndex))% Giant loop to deal with data from each trial
                
                %indicates the overall trial Number
                overallCounter = j;
                
                % Set default for each trial as no dropped samples
                dropThisTrial = 0;
                dropPoint = 0;
                dropPointIndex = 0;
                
                % if the trial doesn't start recording until > 300 ms into the trial, drop this trial (sampling problem with tracker)
                if trialStartTime(j) > .3 
                    dropThisTrial = 1;
                else
                    % Check to make sure no samples were dropped; if they were, check to see if they were during movement
                    
                    % Look for dropped samples.   If two consecutive rows differ by more than 1 in the sample column, this means at least 
                    % 1 sample was dropped.
                    timePointDifferences = diff(movementData(trialStartIndex(j):trialEndIndex(j),sampleMarker));
                    
                    % If more than 5 samples were dropped in trial...
                    if any(timePointDifferences > 10)
                        % As long as that didn't happen twice and the drop was not more than 75 samples (300 ms)...
                        if length(find(timePointDifferences > 10)) < 2 && ~any(timePointDifferences > 300)
                            % Check to see if movement occurred during the drop & find point where the drop occurred
                            dropPoint = find(timePointDifferences > 10);
                            % Calculcate the index # in the file where the drop occurred
                            dropPointIndex = dropPoint + trialStartIndex(j);
                            % Calculate whether the z-index moved more than one in during the drop by calculating the distance moved
                            % between when the drop started and when the drop ended + a few trials (the minimum # samples for threshold)
                            dropPointMovement = diff(movementData(dropPointIndex-1:dropPointIndex+1, z1Marker));
                            
                        else % If two conditions mentioned above weren't met, drop the trial
                            dropThisTrial = 1;
                        end
                    end
                end
                
                % IF it's a dropped trial, mark dropped trial counter, set reach start and end indeces to 1
                if dropThisTrial
                    droppedTrials(overallCounter) = 1;
                    reachStartIndex(j) = 1;
                    reachEndIndex(j) = 1;
                else
                    % mark dropped trial counter
                    droppedTrials(overallCounter) = 0;
                    % see where in the trial min # velocity samples were first exceeded
                    temp = minSamplesExceeded((minSamplesExceeded > trialStartIndex(j) & minSamplesExceeded < trialEndIndex(j)));
                    % If any instances were found, that is where the reach was started
                    if temp
                        reachStartIndex(j) = min(temp);
                        % If no instances were found, there was no movement for this trial; set reach start and end indices to 1
                    else
                        reachStartIndex(j) = 1;
                        reachEndIndex(j) = 1;
                        noMovement(overallCounter) = 1;
                    end
                    
                    % Double check that if no response was recorded, the trial wasn't marked as "correct" in the data file.  If it was,
                    % program will end and return this error message
                    if ~temp
                        if acc(overallCounter)
                            'Error - No velocity threshold exceeded but correct answer recorded!'
                            return;
                        end
                    end
                    
                    % See where in the trial the velocity went back below threshold after the reach started
                    temp2 = minEndSamplesExceeded((minEndSamplesExceeded > reachStartIndex(j) & minEndSamplesExceeded <...
                        trialEndIndex(j)));
                    % If a movement was started...
                    if ~noMovement(overallCounter)
                        % If an end to the movement was discovered...
                        if temp2
                            % Calculate the reach end index as the point where velocity went below threshold for X consecutive samples
                            % (minEndSamples).  The end of the movement is the last of those samples, not the first
                            reachEndIndex(j) = min(temp2);
                            % If no end of the movement was discovered...
                        else
                            % Set the end index as the end of the current trial. Assuming it never went below threshold, this could 
                            % happen if, for example, the movement was still at relatively high speed as the trial ended
                            reachEndIndex(j) = trialEndIndex(j);
                        end
                    end
                end
                
                % If re-using thresholds, figure out if this is a trial that's supposed to be discarded
                if reUseThresholds && ~isempty(dataMatrix2)
                    if dropsCounter <= length(discardBlock)
                        if i == discardBlock(dropsCounter) && j == discardTrial(dropsCounter)
                            'out!'
                            droppedTrials(overallCounter) = 1;
                            dropsCounter = dropsCounter+1;
                            reachStartIndex(j) = 1;
                            reachEndIndex(j) = 1;
                        end
                    end
                end
                            
                % Here start process of calculating curvature
                if ~noMovement(overallCounter)  && ~droppedTrials(overallCounter) %Only bother calculating curvature
                    %if this trial had a movement and if it's not a dropped trial
                    if clickedThreshold
                        reachStartIndex(j) = trialStartIndex(j) + newStartIndex;
                        reachEndIndex(j) = trialStartIndex(j) + newEndIndex ;
                    end
                    
                    % If re-using thresholds, see if this is a trial where we should re-use a threshold
                    if reUseThresholds
                        showTrial = 0;
                        if ~isempty(dataMatrix)
                            if changesCounter <= length(changesBlock)
                                if  (changesBlock(changesCounter) == i && changesTrial(changesCounter) == j)
                                    while (changesCounter < length(changesBlock) & ((changesBlock(changesCounter) == ...
                                            changesBlock(changesCounter + 1) && changesTrial(changesCounter) == ...
                                            changesTrial(changesCounter + 1))))
                                        changesCounter = changesCounter + 1;
                                    end
                                    %'executed'
                                    reachStartIndex(j) = trialStartIndex(j) + newStartIndices(changesCounter);
                                    reachEndIndex(j) = trialStartIndex(j) + newEndIndices(changesCounter);
                                    showTrial = 1;
                                    changesCounter = changesCounter + 1;
                                end
                            end
                        end
                    end
                    
                    % calculate start and end x and y positions on each trial
                    startx1(j) = movementData(reachStartIndex(j),x1Marker);
                    starty1(j) = movementData(reachStartIndex(j),y1Marker);
                    startz1(j) = movementData(reachStartIndex(j),z1Marker);
                    endx1(j) = movementData(reachEndIndex(j),x1Marker);
                    endy1(j) = movementData(reachEndIndex(j),y1Marker);
                    endz1(j) = movementData(reachEndIndex(j),z1Marker);
                    
                    startx2(j) = movementData(reachStartIndex(j),x2Marker);
                    starty2(j) = movementData(reachStartIndex(j),y2Marker);
                    startz2(j) = movementData(reachStartIndex(j),z2Marker);
                    endx2(j) = movementData(reachEndIndex(j),x2Marker);
                    endy2(j) = movementData(reachEndIndex(j),y2Marker);
                    endz2(j) = movementData(reachEndIndex(j),z2Marker);
                    
                    
                    % if you want to avoid looking at the first or last few measurements for curvature because they might be skewed, add 
                    % in a buffer.  Set to zero by default
                    buffer = 0;
                    if reachStartIndex(j) - reachEndIndex(j) <= 10
                        buffer = 0;
                    end
                    indexStartPoint = reachStartIndex(j);
                    indexEndPoint = reachEndIndex(j);
                    
                    % calculate each x and y position throughout the movement
                    currentX1 = movementData(indexStartPoint:indexEndPoint, x1Marker);
                    currentY1 = movementData(indexStartPoint:indexEndPoint, y1Marker);
                    currentZ1 = movementData(indexStartPoint:indexEndPoint, z1Marker);
                    
                    currentX2 = movementData(indexStartPoint:indexEndPoint, x2Marker);
                    currentY2 = movementData(indexStartPoint:indexEndPoint, y2Marker);
                    currentZ2 = movementData(indexStartPoint:indexEndPoint, z2Marker);
                    
                    % Code to calculate the initial angle of deviation up through X% of the movement
                    XPercent = .2;
                    tempSize = size(currentX1,1);
                    indexAtXPercent = round(tempSize * XPercent);
                    maxXPercentPoint = indexStartPoint + indexAtXPercent - 1;
                    
                    tempSize2 = size(currentX2,1);
                    indexAtXPercent2 = round(tempSize2 * XPercent);
                    maxXPercentPoint2 = indexStartPoint + indexAtXPercent2 - 1;
                    
                    currentXPercentPoint = movementData(maxXPercentPoint,x1Marker);
                    currentYPercentPoint = movementData(maxXPercentPoint,y1Marker);
                    
                    currentXPercentPoint2 = movementData(maxXPercentPoint2,x2Marker);
                    currentYPercentPoint2 = movementData(maxXPercentPoint2,y2Marker);
                    
                    % calculate three lines for each position during the movement to form a series of triangles -
                    % start to finish, start to current, and current to finish
                    startToEndLine = pdist2([startx1(j), starty1(j)], [endx1(j), endy1(j)]);
                    startToCurrentLine = pdist2([currentX1, currentY1], [startx1(j), starty1(j)]);
                    currentToEndLine = pdist2([currentX1, currentY1], [endx1(j), endy1(j)]);
                    startToEndLine2 = pdist2([startx2(j), starty2(j)], [endx2(j), endy2(j)]);
                    startToCurrentLine2 = pdist2([currentX2, currentY2], [startx2(j), starty2(j)]);
                    currentToEndLine2 = pdist2([currentX2, currentY2], [endx2(j), endy2(j)]);
                    
                    %Trying to get intial deviation angle, at X% of movement
                    startToXLine =  pdist2([startx1(j), starty1(j)], [currentXPercentPoint, currentYPercentPoint]);
                    XtoEndLine = pdist2([currentXPercentPoint, currentYPercentPoint],[endx1(j), endy1(j)]);
                    startToXLine2 =  pdist2([startx2(j), starty2(j)], [currentXPercentPoint, currentYPercentPoint]);
                    XtoEndLine2 = pdist2([currentXPercentPoint, currentYPercentPoint],[endx2(j), endy2(j)]);
                    
                    % Use law of cosines to get angle of that triangle, which is the angle of deviation from the main movement path 
                    % (defined as a line connecting the start and end points) for each point
                    cosA = (startToCurrentLine.^2 + startToEndLine.^2 - currentToEndLine.^2)./(2*startToCurrentLine*...
                        startToEndLine);
                    angularDeviation = abs(acosd(cosA));
                    
                    cosA2 = (startToCurrentLine2.^2 + startToEndLine2.^2 - currentToEndLine2.^2)./(2*startToCurrentLine2*...
                        startToEndLine2);
                    angularDeviation2 = abs(acosd(cosA2));
                    
                    cosXPercent= (startToXLine^2 + startToEndLine^2 - XtoEndLine^2)/(2*startToXLine*startToEndLine);
                    angularXDeviation(overallCounter) = abs(acosd(cosXPercent));
                    
                    cosXPercent2= (startToXLine2^2 + startToEndLine2^2 - XtoEndLine2^2)/(2*startToXLine2*startToEndLine2);
                    angularXDeviation2(overallCounter) = abs(acosd(cosXPercent2));
                    
                    % calculate the distance (in cm) of deviation from the main movement path at the each point.  Do this by creating 
                    % triangle with angle and current line.  Then use Sin function to calculate opposite line, which would connect from 
                    % current point to startToEndLine
                    distanceFromPath = abs(sind(angularDeviation) .* startToCurrentLine);
                    distanceFromPath2 = abs(sind(angularDeviation2) .* startToCurrentLine2);
                    
                    % Now, use Cos function to calculate adjacent line, which is corollary of startToEndLine.
                    lengthOfAdjacentLine = abs(cosd(angularDeviation) .* startToCurrentLine);
                    lengthOfAdjacentLine2 = abs(cosd(angularDeviation2) .* startToCurrentLine2);
                    
                    % Now, if any of these lengths are longer than the length of the actual startToEndLine, replace distance from path
                    % with distance between current point and end point. This way, if a movement is curved and passes by its target,
                    % we will have a more accurate measure of deviation
                    overShootIndex = find(lengthOfAdjacentLine > startToEndLine);
                    distanceFromPath(overShootIndex) = currentToEndLine(overShootIndex);
                    
                    overShootIndex2 = find(lengthOfAdjacentLine2 > startToEndLine2);
                    distanceFromPath2(overShootIndex2) = currentToEndLine2(overShootIndex2);
                   
                    % Find the highest value for each measure of deviation (angle and distance) for each movement. Divide distance 
                    % measurement by amplitude of movement to normalize
                    maxAngularDeviation(overallCounter) = max(angularDeviation);
                    maxDeviationUncorrected(overallCounter) = max(distanceFromPath);
                    maxDeviation(overallCounter) = max(distanceFromPath)/startToEndLine;
                    
                    maxAngularDeviation2(overallCounter) = max(angularDeviation2);
                    maxDeviationUncorrected2(overallCounter) = max(distanceFromPath2);
                    maxDeviation2(overallCounter) = max(distanceFromPath2)/startToEndLine2;
                    
                    % Find  x and y coordinates for point of maximum curvature
                    maxPoints = find(max(distanceFromPath) == distanceFromPath);
                    maxPoints = maxPoints(1);
                    maxDevPoint(j) = (reachStartIndex(j)+ maxPoints);
                    maxCurrentX(j) = currentX1(maxDevPoint(j) - reachStartIndex(j));
                    maxCurrentY(j) = currentY1(maxDevPoint(j) - reachStartIndex(j));
                    
                    % Calculate index point for sample 1 after max deviation is reached (i.e. when movement correction starts)
                    turnAroundPoint(overallCounter) = maxDevPoint(j) + 1;
                    turnAroundPoint2(overallCounter) = maxDevPoint(j) + 1;
                    
                    % Calculate turnAround RT (i.e. how long into trial did movement re-direct towards target)
                    turnAroundRT(overallCounter) = movementData(turnAroundPoint(overallCounter),timeMarker);
                    turnAroundRT2(overallCounter) = movementData(turnAroundPoint2(overallCounter),timeMarker);
                    
                    % SECTION BELOW HASN'T BEEN USED IN ANY REAL CAPACITY - DOUBLE CHECK BEFORE ASSUMING IT'S RIGHT!
                    
                    %  Calculate deviation at peak velocity
                    peakVelocity(overallCounter) = max(allVelocity(indexStartPoint:indexEndPoint));
                    temp = find(allVelocity(indexStartPoint:indexEndPoint) == peakVelocity(overallCounter));
                    peakVelocityIndex(j) = temp(1);
                    deviationAtPeakVelocity(overallCounter) = distanceFromPath(peakVelocityIndex(j,1))/startToEndLine;
                    
                    %  Calculate deviation at peak acceleration
                    peakAcceleration(overallCounter) = max(allAcceleration(indexStartPoint:indexEndPoint));
                    temp = find(allAcceleration(indexStartPoint:indexEndPoint) == peakAcceleration(overallCounter));
                    peakAccelerationIndex(j) = temp(1);
                    deviationAtPeakAcceleration(overallCounter) = distanceFromPath(peakAccelerationIndex(j,1))/startToEndLine;
                    
                    %  Calculate deviation at peak deceleration
                    peakDeceleration(j) = min(allAcceleration(indexStartPoint:indexEndPoint));
                    temp = find(allAcceleration(indexStartPoint:indexEndPoint) == peakDeceleration(j));
                    peakDecelerationIndex(j) = temp(1);
                    deviationAtPeakDeceleration(overallCounter) = distanceFromPath(peakDecelerationIndex(j,1))/startToEndLine;
                    
                    % Calculate whether curvature was rightward or leftward of main movement line
                    mainVectorAngle = cart2pol((endx1(j) - startx1(j)), (endy1(j) - starty1(j)));
                    curvedVectorAngle = cart2pol((maxCurrentX(j) - startx1(j)), (maxCurrentY(j) - starty1(j)));
                    XVectorAngle = cart2pol((currentXPercentPoint - startx1(j)), (currentYPercentPoint - starty1(j)));
                    
                    if abs(mainVectorAngle) > abs(curvedVectorAngle)
                        rightward(overallCounter) = 1;
                    else
                        rightward(overallCounter) = 0;
                    end
                    
                    if abs(mainVectorAngle) > abs(XVectorAngle)
                        rightwardX(overallCounter) = 1;
                    else
                        rightwardX(overallCounter) = 0;
                    end
                    
                    reachLength = reachEndIndex(j) - reachStartIndex(j) + 1;
                    AllXMovementPoints((overallCounter),1:(reachLength)) = movementData(reachStartIndex(j):reachEndIndex(j),...
                        x1Marker);
                    AllXMovementPoints((overallCounter),reachLength + 1:end) = 0;
                    AllYMovementPoints((overallCounter),1:(reachLength)) = -movementData(reachStartIndex(j):reachEndIndex(j),...
                        y1Marker);
                    AllYMovementPoints((overallCounter),reachLength + 1:end) = 0;
                    
                    % Section to draw images
                    if drawImages &&  (startMiddle == 0 || (j >= trialEnd))
                        
                        clickedThreshold = 0;
                        
                        % Call the figure for movement direction
                        figure(h);
                        hold on
                        clf;
                        
                        % Draw the stimuli on the screen;
                        allData = load([subjDir,filesep,num2str(subjNum),'block',num2str(b),'_grp_all.mat']);
                        eccPx = allData.eccPx;
                        tfDistPx = allData.tfDistPx;
                        barLenPx = allData.barLenPx;
                        hemiIndex = allData.hemiIndex;
                        xCen = allData.xCen;
                        yCen = allData.yCen;
                        
                        clear pos
                        pos = [xCen + eccPx*hemiIndex(j), yCen*-1]; % center bar
                        
                        if targetID(j)
                            pos(2,:) = [xCen + eccPx*hemiIndex(j) + tfDistPx/sqrt(2), (yCen - tfDistPx/sqrt(2))*-1]; % UR
                            pos(3,:) = [xCen + eccPx*hemiIndex(j) - tfDistPx/sqrt(2), (yCen + tfDistPx/sqrt(2))*-1]; % LL
                            pos(4,:) = [xCen + eccPx*hemiIndex(j) - tfDistPx/sqrt(2), (yCen - tfDistPx/sqrt(2))*-1]; % UL
                            pos(5,:) = [xCen + eccPx*hemiIndex(j) + tfDistPx/sqrt(2), (yCen + tfDistPx/sqrt(2))*-1]; % LR
                        end
                        
                        for idx = 1:size(pos,1)
                            hold on
                            plot(pos(idx,1), pos(idx,2), 'bo', 'MarkerSize', 5)
                        end
                        
                        if baseOri(j) % IF HORIZONTAL
                            pos_1 = [pos(1,1)+barLenPx/2, pos(1,2)];
                            pos_2 = [pos(1,1)-barLenPx/2, pos(1,2)];
                        else % IF VERTICAL
                            pos_1 = [pos(1,1), pos(1,2)-barLenPx/2];
                            pos_2 = [pos(1,1), pos(1,2)+barLenPx/2];
                        end
                        plot([pos_1(1),pos_2(1)],[pos_1(2),pos_2(2)],'b-');
                        
                        % Draw the reach movement
                        for n = reachStartIndex(j):reachEndIndex(j)
                            if n == reachEndIndex(j) % Mark the end of the movement with red & green Xs
                                plot (movementData(n,4), -movementData(n,5),  'rx','MarkerSize', 10);
                                plot (movementData(n,6), -movementData(n,7),  'gx','MarkerSize', 10);
                            elseif n == reachStartIndex(j) % Mark the start of the movement with a blue X
                                plot (movementData(n,4), -movementData(n,5),  'bx','MarkerSize', 10);
                                plot (movementData(n,6), -movementData(n,7),  'kx','MarkerSize', 10);
                            elseif n == (maxDevPoint(j)) % Mark the point of max deviation with a black x
                                plot (movementData(n,4), -movementData(n,5),  'kx','MarkerSize', 20);
                                plot (movementData(n,6), -movementData(n,7),  'bx','MarkerSize', 20);
                            else % Mark everything else with a green O
                                if reachedTo((overallCounter))==4 || reachedTo((overallCounter))==1
                                    plot (movementData(n,4), -movementData(n,5),'go');
                                    plot (movementData(n,6), -movementData(n,7),'ro');
                                elseif  reachedTo((overallCounter))==2
                                    plot (movementData(n,4), -movementData(n,5),'go');
                                end
                            end
                            hold on;
                        end
                        
                        % Set the axis so its consistent for each trial
                        axis([0 1280 -1200 0]);
                        figure(h);
                        tempTitle = (strcat('Subject#', num2str(subjNum),' Trial#',num2str(j), ' Block#', num2str(b), ...
                            ', MaxCurvature=', num2str(maxDeviation(overallCounter) *signOfCurvature(overallCounter)),...
                            ', Key Acc=', num2str(acc(overallCounter)), ', Grasp Acc=',num2str(acc_grasp(overallCounter))));
                        title(tempTitle);
                        hold off;
                        
                        % Create a velocity profile for each trial
                        velocityProfile = allVelocity(trialStartIndex(j):trialEndIndex(j));
                        trialsample(j) = length(trialStartIndex(j):trialEndIndex(j));
                        trialtime (j) = movementData(trialEndIndex(j),timeMarker)-movementData(trialStartIndex(j),timeMarker);
                        figure(hh);
                        ha = area([reachStartIndex(j) - trialStartIndex(j),reachEndIndex(j) - trialStartIndex(j)],[2000,...
                            2000]);
                        hold on;
                        t(j)=movementData(reachStartIndex(j),timeMarker);
                        tempTitle2 = (strcat( 'starttime:', num2str(t(j))));
                        title(tempTitle2);
                        plot(velocityProfile); % orange plot
                        axis([0 700 0 150]);
                        figure(hh);
                        drawnow;
                        breakloop = 0;
                        loopStart = tic;
                        while ~breakloop && (toc(loopStart) < timeLimit)
                            [keyIsDown, secs, keyCode]= KbCheck;
                            key = KbName(keyCode);
                            if keyIsDown
                                keyPressed = key(1);
                                switch keyPressed
                                    case 'b' % go back
                                        if j == previousj
                                            j = previousj - 4;
                                        else
                                            j = previousj - 1;
                                        end
                                        if j < 0
                                            j = 1;
                                        end
                                        breakloop = 1;
                                        
                                    case 'd'
                                        ListenChar(0);
                                        prompt{1} = 'Notes';
                                        default{1} = 'Dropped Samples';
                                        lines = ones(1,1);
                                        answer = inputdlg(prompt, 'Notes', lines, default);
                                        note = answer{1};
                                        prompt2{1} = 'Undo?';
                                        default2{1} = '1';
                                        lines2 = ones(1,1);
                                        answer2 = inputdlg(prompt2, 'Undo?', lines2, default2);
                                        undoResp = str2num(answer2{1});
                                        if undoResp == 0
                                            droppedTrials(overallCounter) = 1;
                                            reachStartIndex(j) = 1;
                                            reachEndIndex(j) = 1;
                                            fprintf(dataFilepointer, '%g\t %g \t  %g \t %s \n', allFiles,i,j,note);
                                        else
                                            j = j-1;
                                        end
                                        ListenChar(2);
                                        breakloop = 1;
                                        
                                    case 't'
                                        clickedThreshold = 1;
                                        ListenChar(0);
                                        [axisX,axisY] = ginput(2);
                                        breakLoop2 = 0;
                                        while ~breakLoop2
                                            [keyIsDown, secs, keyCode]= KbCheck;
                                            key = KbName(keyCode);
                                            if keyIsDown
                                                keyPressed = key(1);
                                                switch keyPressed
                                                    case '/'
                                                        newStartIndex = round(axisX(1));
                                                        newEndIndex = round(axisX(2));
                                                        breakLoop2 = 1;
                                                    case '1'
                                                        newStartIndex = round(axisX(1));
                                                        newEndIndex = reachEndIndex(j) - trialStartIndex(j);
                                                        breakLoop2 = 1;
                                                    case '9'
                                                        newStartIndex = reachStartIndex(j) - trialStartIndex(j);
                                                        newEndIndex = round(axisX(2));
                                                        breakLoop2 = 1;
                                                    otherwise
                                                        WaitSecs(.1);
                                                end
                                            end
                                        end
                                        
                                        fprintf(dataFilepointer2, '%g\t %g \t  %g \t  %g \t %g \n', allFiles,i,j,...
                                            newStartIndex, newEndIndex);
                                        
                                        j = j - 1;
                                        
                                        ListenChar(2);
                                        breakloop = 1;
                                        
                                    case 'z'
                                        breakloop = 1;
                                        clickedThreshold = 0;
                                        
                                    case 's'
                                        timeLimit = 10000;
                                    case 'f'
                                        timeLimit = originalTimeLimit;
                                        clickedThreshold = 0;
                                    case 'j'
                                        i = i - 1;
                                        j = 0;
                                        breakloop = 1;
                                        
                                    case 'q'
                                        cd ..
                                        cd ..
                                        close all;
                                        ListenChar(0);
                                        return;
                                    otherwise
                                        breakloop = 0;
                                end
                            end
                        end
                        previousj = j;
                        clf;
                    end
                    
                    
                    
                end
                
                % increment the loop
                j = j+1;
            end
            
            % Calculate Response Time (how long did it take to initiate a response) and 
            % Movement Time (how long did it take to reach the target) for each trial
            movementStartTime = movementData(reachStartIndex, 3);
            movementEndTime = movementData(reachEndIndex, 3);
            %         RT(1+((i-1) * blockLength):blockLength+((i-1) * blockLength)) = movementStartTime - trialStartTime;
            %         MT(1+((i-1) * blockLength):blockLength+((i-1) * blockLength))  = movementEndTime - movementStartTime;
            RT = movementStartTime - trialStartTime;
            MT = movementEndTime - movementStartTime;
            
        end
        
        cd(subjDir)
        close all
        if saveFiles
            backupFileName = ['MovementData_', curSubj, '.mat'];
            backupFile = [subjDir, filesep, backupFileName];
            save(backupFile);
        end
        
    end
    
    ListenChar(0);
    close all
end

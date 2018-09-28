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

%subject = input('\n Sub number: ');
%subSeq = input('\n Which number of new subject is this?  ');
% subject=[];   % subject number
% Totalsub = 13;      
% subSeq=[];    % saved file name for each subject
 subject=[12,21,22,23,24,32,33,35,36,37,38,39,40,41,43,44];   % subject number
 Totalsub = 16;      
 subSeq=[12,21,22,23,24,32,33,35,36,37,38,39,40,41,43,44];    % saved file name for each subject
% Instructions.  If examining images and adjusting thresholds by hand, use
% the following keys:
% 's' - stay at this image
% 'f' - foward to following images
% 'z' - skip to next image, current is OK
% 'b' - go back one image
% 'd' - delete the current image, leave a note as to why.  will be saved to file
% 't' - adjust threshold by hand to see if it makes the movement make more
% sense.  This will be written to file.  Please note, if you do multiple
% threshold adjustments for the same trial, you have to manually delete the
% earlier ones in the file afterwards (unless I fix this).
% you have to click "return" after clicking the new beginning and end
% thresholds.  Then, click "return" again to confirm, or click "1" if you
% want to only change the beginning threshold (it will ignore where you
% clicked the end threshold), or "9" if you only want to change the end
% threshold.

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
% allVars(1).name = 'TrialNumber';
% allVars(1).number = '1';
% allVars(2).name = 'DistractorPresence';
% allVars(2).number = '5';
%[ reachedTo,timeElapsed, acc(i), targetLocation(i), targetColor(i),  distractorLocation(i), distractorColor(i),targetShape(i), trialType(i)],'-append', 'roffset', [],'delimiter', '\t');
allVars(1).name = 'targetLocation';
allVars(1).number = '3';
allVars(2).name = 'r1';
allVars(2).number = '4';
allVars(3).name = 'r2';
allVars(3).number = '5';
allVars(4).name = 'reachto';
allVars(4).number = '6';
allVars(5).name = 'timeElapsed';
allVars(5).number = '9';
allVars(6).name = 'ikeyacc';
allVars(6).number = '11';
allVars(7).name = 'iacc';
allVars(7).number = '12';
allVars(8).name = 'markerWait';
allVars(8).number = '13';

orip=pwd;

for numsub = 1:Totalsub
    numsub
    cd(orip);
    currentSubjectNumber = subject(numsub);%Number of subject you are currently analyzing.  
    % Note that this just assigns numbers to the subject folders in
    % alphabetical order, so if you're first subject is #305, that will be
    % subject #1 for these purposes.
    % Various flags
    drawImages = 1; %1 if you want to see each threshold
    saveFiles = 0;%1 if you want to save the new files
    reUseThresholds = 1;%1 if you want to re-use thresholds from a previous pre-processing
    startMiddle = 0;% 1 if you want to start in the middle - this is uesful if you crash
    %in the middle of pre-processing, so you don't have to start from the
    %beginning.  It will prompt you asking which block and trial number you
    %want to start at
    drawInPixels = 1;% currently only works as 1
    originalTimeLimit = .2; % How long should each image be displayed before auto-advancing
    timeLimit = originalTimeLimit;
    % Specify the length of each block in the current study
    blockLength = 32;
    numBlocks = 12;
    blockLength = blockLength * numBlocks;
    % mark as "1" if you want to upsample data
    upSample = 0;
    samplingRate = 200;

    % Current subject for analysis, if running only one subject at a time.
    for currentSub = currentSubjectNumber
        %currentSub = 14;

        if startMiddle
            blockEnd = str2num(input('\n Enter last Block Number: ', 's'));
            trialEnd = str2num(input('\n Enter last Trial Number: ', 's'));
        end           


        % Name of folder where you'll keep change files
        if ~isdir('ManualChanges')
            mkdir('ManualChanges');
        end
        subFolder = 'ManualChanges/exp1';

        % Name of folder where data is stored
        dataFolder = '../Subject_Folder/exp1';

        % Name of processed .mat files to be created
        matFileName = 'MovementData_';

%         if ~isdir('ProcessedFiles')
%             mkdir('ProcessedFiles');
%         end
        % Where those processed .mat files should be stored
        matFilePath = 'Subject_Folder/exp1';

        ListenChar(2);%suppress keyboard output in the command line

        if saveFiles && ~reUseThresholds
            % % Create name for output file for Redo Trials (trials where visual
            % inspection suggests something is wrong) and open it
            outputFilename = strcat(subFolder,filesep,'RedoTrials',num2str(currentSub));
            dataFilename = strcat(outputFilename, '.txt');
            dataFilepointer = fopen(dataFilename, 'w');
            %fprintf(dataFilepointer, '%s \t  %s \t %s \t  %s \n','Subject','Block', 'Trial', 'Notes');

            % Create name for output file for modified threshold Trials (trials where the movement thresholds were manually changed) and open it
            outputFilename2 = strcat(subFolder,filesep,'ThresholdAdjustments',num2str(currentSub));
            dataFilename2 = strcat(outputFilename2, '.txt');
            dataFilepointer2 = fopen(dataFilename2, 'w');
            %fprintf(dataFilepointer2, '%s \t  %s \t %s \t %s \t %s \n','Subject','Block', 'Trial', 'NewStartIndex', 'NewEndIndex');
        elseif reUseThresholds %&& startMiddle
            % % Create name for output file for Redo Trials (trials where visual
            % inspection suggests something is wrong) and open it
            outputFilename = strcat(subFolder,filesep,'RedoTrials',num2str(currentSub));
            dataFilename = strcat(outputFilename, '.txt');
            dataFilepointer = fopen(dataFilename, 'a');
            %fprintf(dataFilepointer, '%s \t  %s \t %s \t  %s \n','Subject','Block', 'Trial', 'Notes');

            % Create name for output file for modified threshold Trials (trials where the movement thresholds were manually changed) and open it
            outputFilename2 = strcat(subFolder,filesep,'ThresholdAdjustments',num2str(currentSub));
            dataFilename2 = strcat(outputFilename2, '.txt');
            dataFilepointer2 = fopen(dataFilename2, 'a');
            %fprintf(dataFilepointer2, '%s \t  %s \t %s \t %s \t %s \n','Subject','Block', 'Trial', 'NewStartIndex', 'NewEndIndex');
        end


        %Change to path of data
        oripath=pwd;
        path= strcat(pwd,filesep,dataFolder,filesep);
        cd (path);
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
        for allFiles = currentSub % Can be changed to pre-process multiple subjects at once.  Should
            %   use only if re-doing analysis with pre-processed files

            % Code to call up files with info about which manual changes have been
            % made
            if reUseThresholds % && startMiddle
                %cd ..  
                changeFile = strcat(oripath,filesep,subFolder,filesep,'ThresholdAdjustments',num2str(currentSub),'.txt');
                discardFile = strcat(oripath,filesep,subFolder,filesep,'RedoTrials',num2str(currentSub),'.txt');
                dataMatrix = load(changeFile);
                dataMatrix2 = load(discardFile);
                if ~isempty(dataMatrix)
                    changesBlock = dataMatrix(:,2);
                    changesTrial = dataMatrix(:,3);
                    newStartIndeces = dataMatrix(:,4);
                    newEndIndeces = dataMatrix(:,5);
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
            curDirectory = strcat(path,'sub_',num2str(subject(numsub)),'/');
            cd (curDirectory);
            files = dir(curDirectory);
            files(1:2) = [];
            fileIndexStart = 1;
            if strcmp(files(1).name,'.DS_Store')
                files(1) = [];
            end                

            % Which column in the data file is used to measure x & y coordinates
            % Currently in space (and thus inches - need to convert to cm)        
            x1Marker = 8;
            z1Marker = 10;
            y1Marker = 9;
            x2Marker = 11;
            z2Marker = 13;
            y2Marker = 12;
            timeMarker = 3;
            trialNumberMarker = 2;
            sampleMarker = 14;

            % This is where the threshold is set for what determines the start
            % and end of a movement.  The min samples portions determine how
            % many consectuive samples must exceed (or fall below) the
            % threshold to determine the start/end of a movement.
            standardThreshold = 15; % cm/s
            standardLowerThreshold = 15; 
            standardMinSamples = 1;
            standardMinEndSamples = 1;


            % Set minimum thresholds (in cm/s) to define the start and end of arm movement,
            % i.e. what velocity (or acceleration, if needed) must be exceeded to consider
            % something an arm movement
            threshold = standardThreshold;
            lowerThreshold = standardLowerThreshold;

            % Set minimun number of consecutive samples over/under threshold
            % required to accept the start/end of a movement
            minSamples = standardMinSamples;
            minEndSamples = standardMinEndSamples;

            % Calculate total number of runs for each subject
            numRuns = 1; % length(files)/2;

            % preallocate variables that will be recorded as totals for all blocks
            totalTrials = blockLength*numRuns;
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

            % for loop to analyze individual data files for current subject
            for i = 1:numRuns                                  

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

                allVelocity = [];

                % Open two files for each run through the loop  - simple data file = file,
                % and movement data file = file2        
                data=[];
                for iii=1:numBlocks
                    file=strcat(curDirectory,'/',num2str(subject(numsub)),'_1_graspExp_block',num2str(iii),'.txt');
                    data2 = load(file);
                    data=[data;data2];
                end
                
                % delete crashed trials
                index = find(diff(data(:,2))==0);
                block = data(index,1);
                trial = data(index,2);
                data(index,:)=[];
                
                file2=strcat(curDirectory,'/',num2str(subject(numsub)),'_1_graspExp.txt');

                if upSample
                    tempmovementData = load(file2);
                    tempa = [1;diff(tempmovementData(:,2))];
                    temptrialNumbers = tempmovementData(:,2);
                    temptrialStartIndex = find(tempa);
                    temptrialEndIndex = [temptrialStartIndex(2:blockLength) - 1;length(tempmovementData)];            
                    movementData = [];            

                    for j = 1:length(temptrialStartIndex)
                       timeElapsed = tempmovementData(temptrialEndIndex(j),timeMarker) - tempmovementData(temptrialStartIndex(j),timeMarker);
                       upsampleRate = round(timeElapsed * samplingRate);
                       numSamples = temptrialEndIndex(j) - temptrialStartIndex(j);               
                       if samplingRate > (numSamples/timeElapsed)                     
                             tempData = zeros(upsampleRate,9);
                             tempData(:,1) = tempmovementData(temptrialStartIndex(j),1);
                             tempData(:,2) = tempmovementData(temptrialStartIndex(j),2);
                             newpath=pwd;
                             cd(oripath);
                             tempData(:,3) = EZResample(tempmovementData(temptrialStartIndex(j):temptrialEndIndex(j),3),upsampleRate,0);
                             tempData(:,4) = EZResample(tempmovementData(temptrialStartIndex(j):temptrialEndIndex(j),4),upsampleRate,0);
                             tempData(:,5) = EZResample(tempmovementData(temptrialStartIndex(j):temptrialEndIndex(j),5),upsampleRate,0);
                             tempData(:,6) = EZResample(tempmovementData(temptrialStartIndex(j):temptrialEndIndex(j),6),upsampleRate,0);
                             tempData(:,7) = EZResample(tempmovementData(temptrialStartIndex(j):temptrialEndIndex(j),7),upsampleRate,0);
                             tempData(:,8) = EZResample(tempmovementData(temptrialStartIndex(j):temptrialEndIndex(j),8),upsampleRate,0);
                             tempData(:,9) = EZResample(tempmovementData(temptrialStartIndex(j):temptrialEndIndex(j),9),upsampleRate,0);
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
                delete_index=[];
                ii = 1;
                for jj = 1:length(block);
                    time=0;
                    flag=0;
                    while flag==0
                        if any(movementData(ii,1)==block(jj))&& any(movementData(ii,2)==trial(jj))&& (movementData(ii,3)>time)
                            delete_index=[delete_index,ii];
                            time = movementData(ii,3);                    
                        elseif any(movementData(ii,1)==block(jj))&& any(movementData(ii,2)==trial(jj))&& (movementData(ii,3)<time)
                            flag=1;
                        end;
                        ii=ii+1;
                    end
                end;
                movementData(delete_index,:)=[];
                
                for j = 1:length(allVars)
                    eval([allVars(j).name,'(((i-1) * blockLength + 1):(i*blockLength)) = data(:,',allVars(j).number,');']);
                end           
                % Calculate speed for every sample in the file. A zero is added at the beginning because point 1 cm the file will not
                % have a speed, and the diff command will output the first speed
                % value, from point 2, as the first value in the scalar.  Adding a
                % zero will line up the vector correctly
                timeChange = diff(movementData(:,timeMarker));        
                xArmChange = -diff(movementData(:,x1Marker)); % Distance between each sample in X dimension
                yArmChange = -diff(movementData(:,y1Marker)); % Distance between each sample in Y dimension
                zArmChange = -diff(movementData(:,z1Marker)); % Distance between each sample in Z dimension
                totalArmChange = sqrt(xArmChange.^2 + yArmChange.^2 + zArmChange.^2); % Calculate the distance (currently 2D) between each movement
                tempAllVelocity = [0;totalArmChange./timeChange]; % Calculate velocity by dividing distance by time change                                  


                % Calculate acceleration for every sample in the file.  Add a zero at the
                % beginning for same logic above - first data point will not have
                % an acceleration
                velocityChange = diff(tempAllVelocity); % Change in velocity between each sample
                allAcceleration = [0;velocityChange./timeChange]; % Calculate acceleration by dividing velocity change by time change

                a = [1;diff(movementData(:,trialNumberMarker))];
                trialNumbers = movementData(:,trialNumberMarker);
                %trialStartIndex = find(a);

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
                        timeElapsedForThisTrial = movementData(trialEndIndex(xx),timeMarker) - movementData(trialStartIndex(xx), timeMarker);
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
                            [b,a] = butter(nthOrder,butterWindow);
                            tempx = filtfilt(b,a,tempAllVelocity(trialStartIndex(xx):trialEndIndex(xx)));
                            allVelocity(trialStartIndex(xx):trialEndIndex(xx)) = tempx;
                        end
                        % The line below is just for use when butterworth
                        % filter is unavailable
                        %allVelocity(trialStartIndex(xx):trialEndIndex(xx)) = tempAllVelocity(trialStartIndex(xx):trialEndIndex(xx));
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
                    overallCounter = j+((i-1) * blockLength);

                    % Set default for each trial as no dropped samples
                    dropThisTrial = 0;
                    dropPoint = 0;
                    dropPointIndex = 0;

                    % if the trial doesn't start recording until > 300 ms into the
                    % trial, drop this trial (sampling problem with tracker)
                    if trialStartTime(j) > .3
                        dropThisTrial = 1;
                    else
                    

                        % Check to make sure no samples were dropped; if they
                        % were, check to see if they were during movement


                        % Look for dropped samples.   If two consecutive rows differ by
                        % more than 1 in the sample column, this means at least 1 sample was
                        % dropped.
                        timePointDifferences = diff(movementData(trialStartIndex(j):trialEndIndex(j),sampleMarker));

                        % If more than 5 samples were dropped in trial...
                        if any(timePointDifferences > 10)
                            % As long as that didn't happen twice and the drop was not
                            % more than 75 samples (300 ms)...
                            if length(find(timePointDifferences > 10)) < 2 && ~any(timePointDifferences > 300)
                                % Check to see if movement occurred during the drop
                                % find point where the drop occurred
                                dropPoint = find(timePointDifferences > 10);
                                % Calculcate the index # in the file where the drop
                                % occurred
                                dropPointIndex = dropPoint + trialStartIndex(j);
                                % Calculate whether the z-index moved more than
                                % one in
                                % during the drop by calculating the distance moved
                                % between when the drop started and when the drop
                                % ended + a few trials (the minimum # samples for
                                % threshold)
                                dropPointMovement = diff(movementData(dropPointIndex-1:dropPointIndex+1, z1Marker));

                                % If there was at least 1 cm of movement, drop the
                                % trial
%                                 if any(dropPointMovement < -1)
%                                     dropThisTrial = 1;
%                                 end

                            else % If two conditions mentioned above weren't met, drop the trial
                                dropThisTrial = 1;
                            end
                        end
                    end

                    % IF it's a dropped trial, mark dropped trial counter, set
                    % reach start and end indeces to 1
                    if dropThisTrial
                        droppedTrials(overallCounter) = 1;
                        reachStartIndex(j) = 1;
                        reachEndIndex(j) = 1;
                    else
                        % mark dropped trial counter
                        droppedTrials(overallCounter) = 0;
                        % see where in the trial min # velocity samples were first
                        % exceeded
                        temp = minSamplesExceeded((minSamplesExceeded > trialStartIndex(j) & minSamplesExceeded < trialEndIndex(j)));
                        % If any instances were found, that is where the reach was
                        % started
                        if temp
                            reachStartIndex(j) = min(temp);
                            % If no instances were found, there was no Movement for this
                            % trial; set reach start and end indeces to 1
                        else
                            reachStartIndex(j) = 1;%trialStartIndex(j);
                            reachEndIndex(j) = 1;
                            noMovement(overallCounter) = 1;
                        end

                        % Double check that if no response was recorded, the trial
                        % wasn't marked as "correct" in the data file.  If it was,
                        % program will end and return this error message
                        if ~temp
                            if acc(overallCounter)
                                'Error - No velocity threshold exceeded but correct answer recorded!'
                                return;
                            end
                        end

                        % See where in the trial the velocity went back below
                        % threshold after the reach started
                        temp2 = minEndSamplesExceeded((minEndSamplesExceeded > reachStartIndex(j) & minEndSamplesExceeded < trialEndIndex(j)));
                        % If a movement was started...
                        if ~noMovement(overallCounter)
                            % If an end to the movement was discovered...
                            if temp2
                                % Calculate the reach end index as the point where velocity went
                                % below threshold for X consecutive samples
                                % (minEndSamples).  The end of the movement is the
                                % last of those samples, not the first
                                reachEndIndex(j) = min(temp2);
                                % If no end of the movement was discovered...
                            else
                                % Set the end index as the end of the current trial.
                                %  Assuming it never went below threshold, This could happen if,
                                % for example, the movement was still at relatively
                                % high speed as the trial ended
                                reachEndIndex(j) = trialEndIndex(j);
                            end
                        end
                    end

                    % If re-using thresholds, figure out if this is a trial
                    % that's supposed to be discarded
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
                        %if this trial had a movement and if it's not a dropped
                        %trial

                        if clickedThreshold
                            reachStartIndex(j) = trialStartIndex(j) + newStartIndex;
                            reachEndIndex(j) = trialStartIndex(j) + newEndIndex ;
                        end

                        % If re-using thresholds, see if this is a trial where
                        % we should re-use a threshold
                        if reUseThresholds
                            showTrial = 0;
                            if ~isempty(dataMatrix)
                                if changesCounter <= length(changesBlock)
                                    if  (changesBlock(changesCounter) == i && changesTrial(changesCounter) == j)
                                        while (changesCounter < length(changesBlock) & ((changesBlock(changesCounter) == changesBlock(changesCounter + 1) && changesTrial(changesCounter) == changesTrial(changesCounter + 1))))
                                            changesCounter = changesCounter + 1;
                                        end
                                        %'executed'
                                        reachStartIndex(j) = trialStartIndex(j) + newStartIndeces(changesCounter);
                                        reachEndIndex(j) = trialStartIndex(j) + newEndIndeces(changesCounter);
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


                        % if you want to avoid looking at the first or last few measurements for curvature because
                        % they might be skewed, add in a buffer.  Set to zero by
                        % default
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
                        % calculate three lines for each position during the movement 
                        % to form a series of triangles -
                        % start to finish, start to current, and current to finish
                        startToEndLine = pdist2([startx1(j), starty1(j)], [endx1(j), endy1(j)]);
                        startToCurrentLine = pdist2([currentX1, currentY1], [startx1(j), starty1(j)]);
                        currentToEndLine = pdist2([currentX1, currentY1], [endx1(j), endy1(j)]);  
                        startToEndLine2 = pdist2([startx2(j), starty2(j)], [endx2(j), endy2(j)]);
                        startToCurrentLine2 = pdist2([currentX2, currentY2], [startx2(j), starty2(j)]);
                        currentToEndLine2 = pdist2([currentX2, currentY2], [endx2(j), endy2(j)]);                    

                        %Trying to get intial deviation angle, at X% of
                        %movement
                        startToXLine =  pdist2([startx1(j), starty1(j)], [currentXPercentPoint, currentYPercentPoint]);
                        XtoEndLine = pdist2([currentXPercentPoint, currentYPercentPoint],[endx1(j), endy1(j)]);
                        startToXLine2 =  pdist2([startx2(j), starty2(j)], [currentXPercentPoint, currentYPercentPoint]);
                        XtoEndLine2 = pdist2([currentXPercentPoint, currentYPercentPoint],[endx2(j), endy2(j)]);                   
                        % Use law of cosines to get angle of that triangle, which is
                        % the angle of deviation from the main movement path (defined as 
                        % a line connecting the start and end points) for each point           
                        cosA = (startToCurrentLine.^2 + startToEndLine.^2 - currentToEndLine.^2)./(2*startToCurrentLine*startToEndLine);
                        angularDeviation = abs(acosd(cosA));

                        cosA2 = (startToCurrentLine2.^2 + startToEndLine2.^2 - currentToEndLine2.^2)./(2*startToCurrentLine2*startToEndLine2);
                        angularDeviation2 = abs(acosd(cosA2));

                        cosXPercent= (startToXLine^2 + startToEndLine^2 - XtoEndLine^2)/(2*startToXLine*startToEndLine);
                        angularXDeviation(overallCounter) = abs(acosd(cosXPercent));

                        cosXPercent2= (startToXLine2^2 + startToEndLine2^2 - XtoEndLine2^2)/(2*startToXLine2*startToEndLine2);
                        angularXDeviation2(overallCounter) = abs(acosd(cosXPercent2));
                        % calculate the distance (in cm) of deviation from the 
                        % main movement path at the each point.  Do this by
                        % creating triangle with angle and current line.  Then
                        % use Sin function to calculate opposite line, which would
                        % connect from current point to startToEndLine
                        distanceFromPath = abs(sind(angularDeviation) .* startToCurrentLine);
                        distanceFromPath2 = abs(sind(angularDeviation2) .* startToCurrentLine2);
                        % Now, use Cos function to calculate adjacent line, which
                        % is corrolary of startToEndLine.  
                        lengthOfAdjacentLine = abs(cosd(angularDeviation) .* startToCurrentLine);
                        lengthOfAdjacentLine2 = abs(cosd(angularDeviation2) .* startToCurrentLine2);
                        % Now, if any of these lengths are longer than the length
                        % of the actual startToEndLine, replace distance from path
                        % with distance between current point and end point.  This
                        % way, if a movement is curved and passes by its target,
                        % we will have a more accurate measure of deviation
                        overShootIndex = find(lengthOfAdjacentLine > startToEndLine);
                        distanceFromPath(overShootIndex) = currentToEndLine(overShootIndex);

                        overShootIndex2 = find(lengthOfAdjacentLine2 > startToEndLine2);
                        distanceFromPath2(overShootIndex2) = currentToEndLine2(overShootIndex2);
                        % Find the highest value for each measure of deviation (angle
                        % and distance) for each movement.  Divide distance measurement
                        % by amplitude of movement to normalize
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

    %                     maxPoints2 = find(max(distanceFromPath2) == distanceFromPath2); 
    %                     maxPoints2 = maxPoints2(1);
    %                     maxDevPoint2(j) = (reachStartIndex(j)+ maxPoints2);
    %                     maxCurrentX2(j) = currentX1(maxDevPoint2(j) - reachStartIndex(j));
    %                     maxCurrentY2(j) = currentY1(maxDevPoint2(j) - reachStartIndex(j));
                        % Calculate index point for sample 1 after max deviation 
                        % is reached (i.e. when movement correction starts) 
                        turnAroundPoint(overallCounter) = maxDevPoint(j) + 1;
                        turnAroundPoint2(overallCounter) = maxDevPoint(j) + 1;
                        % Calculate turnAround RT (i.e. how long into trial did
                        % movement re-direct towards target)
                        turnAroundRT(overallCounter) = movementData(turnAroundPoint(overallCounter),timeMarker);
                        turnAroundRT2(overallCounter) = movementData(turnAroundPoint2(overallCounter),timeMarker);
                        %COMPointInMovement(overallCounter) = maxPoints;
                        %COMpercentageOfX(overallCounter) = abs(maxCurrentX(j) - startx(j))/abs(endx(j) - startx(j));
                        %COMpercentageOfY(overallCounter) = abs(maxCurrentY(j) - starty(j))/abs(endy(j) - starty(j));

                        % SECTION BELOW HASN'T BEEN USED IN ANY REAL CAPACITY -
                        % DOUBLE CHECK BEFORE ASSUMING IT'S RIGHT!

                        %  Calculate deviation at peak velocity
                        peakVelocity(overallCounter) = max(allVelocity(indexStartPoint:indexEndPoint));
                        temp = find(allVelocity(indexStartPoint:indexEndPoint) == peakVelocity(overallCounter));
                        peakVelocityIndex(j) = temp(1);
                        deviationAtPeakVelocity(overallCounter) = distanceFromPath(peakVelocityIndex(j,1))/startToEndLine;
                        %angleAtPeakVelocity(overallCounter) = angularDeviation(peakVelocityIndex(j,1));

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

                        % Calculate whether curvature was rightward or leftward of 
                        % main movement line                
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
                        AllXMovementPoints((overallCounter),1:(reachLength)) = movementData(reachStartIndex(j):reachEndIndex(j), x1Marker);
                        AllXMovementPoints((overallCounter),reachLength + 1:end) = 0;
                        AllYMovementPoints((overallCounter),1:(reachLength)) = -movementData(reachStartIndex(j):reachEndIndex(j), y1Marker);
                        AllYMovementPoints((overallCounter),reachLength + 1:end) = 0;



                        % Section to draw images
                        if drawImages &&  (startMiddle == 0 || ( j>=trialEnd )) %  &&  ~droppedTrials(overallCounter) && acc(overallCounter) &&  (startMiddle == 0 || ((i > blockEnd || (i == blockEnd && j > trialEnd - 5))))                        
                               
                               clickedThreshold = 0;

                               % Call the figure for movement direction
                               figure(h);
                               clf;

                               % Draw the stimuli on the screen;

                                % SPECIAL SECTION
                               % distance of each object from fixation
                               distFromCenter = 237.22;
                               xCen = 640;
                               yCen = 512;
                               vsize=56;
                               pos(1,:) = [xCen-distFromCenter,yCen];%7:30
                               pos(2,:) = [xCen+distFromCenter,yCen];% 10:30
                               %pos(3,:) = [xCen+distFromCenter,yCen-distFromCenter];% 1:30
                               %pos(4,:) = [xCen+distFromCenter,yCen+distFromCenter];% 4:30
                               rect_cen_x = [pos(1,1), pos(2,1)];
                               rect_cen_y = [pos(1,2), pos(2,2)];

                               for k = 1:2
                                   if k == targetLocation(j+((i-1) * blockLength))
                                       plot(pos(k,1), -pos(k,2), 'bo', 'MarkerSize', vsize./sqrt(2));
                                       pp=k;
                                    %else
                                    %    plot(pos(k,1), -pos(k,2), 'gd', 'MarkerSize', 30);
                                    end
                                   hold on;
                               end

                               if  r2(j+((i-1) * blockLength))==45
                                   pos_top = [pos(pp,1)-vsize./sqrt(2), pos(pp,2)-vsize./sqrt(2)];
                                   pos_bottom = [pos(pp,1)+vsize./sqrt(2), pos(pp,2)+vsize./sqrt(2)];
                                   plot([pos_top(1),pos_bottom(1)],[-pos_top(2),-pos_bottom(2)],'b-');
                               else
                                   pos_top = [pos(pp,1)+vsize./sqrt(2), pos(pp,2)-vsize./sqrt(2)];
                                   pos_bottom = [pos(pp,1)-vsize./sqrt(2), pos(pp,2)+vsize./sqrt(2)];
                                   plot([pos_bottom(1),pos_top(1)],[-pos_bottom(2),-pos_top(2)],'b-');
                               end
                               hold on;
                                % Draw the reach movement 
                               for n = reachStartIndex(j):reachEndIndex(j)                             
                                    if n == reachEndIndex(j) % Mark the end of the movement with a red X
                                        plot (movementData(n,4), -movementData(n,5),  'rx','MarkerSize', 10);  
                                        plot (movementData(n,6), -movementData(n,7),  'gx','MarkerSize', 10);
                                    elseif n == reachStartIndex(j) % Mark the start of the movement with a blue X
                                        plot (movementData(n,4), -movementData(n,5),  'bx','MarkerSize', 10);
                                        plot (movementData(n,6), -movementData(n,7),  'kx','MarkerSize', 10);
                                    elseif n == (maxDevPoint(j)) % Mark the point of max deviation with a black x
                                        plot (movementData(n,4), -movementData(n,5),  'kx','MarkerSize', 20);
                                        plot (movementData(n,6), -movementData(n,7),  'bx','MarkerSize', 20);
                                    else % Mark everything else with a green O
                                        if reachto((overallCounter))==4 || reachto((overallCounter))==1
                                            plot (movementData(n,4), -movementData(n,5),'go');
                                            plot (movementData(n,6), -movementData(n,7),'ro');
                                        elseif  reachto((overallCounter))==2
                                            plot (movementData(n,4), -movementData(n,5),'go');
                                        end
                                    end
                                    hold on;                                                           
                               end

                               % Set the axis so its consistent for each
                               % trial
                               axis([0 1280 -1200 0]);
                               figure(h);
                                % Create title for the movement chart
                               tempTitle = (strcat('Subject #', num2str(allFiles),'Trial #',num2str(j), '__', num2str(i), '___', 'MaxCurvature', num2str(maxDeviation(overallCounter) *signOfCurvature(overallCounter)), '___', 'accuracy', num2str(acc(overallCounter))));
                               title(tempTitle);
                               hold off;
                               % Create a velocity profile for each trial
                               velocityProfile = allVelocity(trialStartIndex(j):trialEndIndex(j));
                               trialsample(j) = length(trialStartIndex(j):trialEndIndex(j));
                               trialtime (j) = movementData(trialEndIndex(j),timeMarker)-movementData(trialStartIndex(j),timeMarker);
                               figure(hh); 
                               ha = area([reachStartIndex(j) - trialStartIndex(j),reachEndIndex(j) - trialStartIndex(j)],[2000, 2000]);
                                hold on;
                                t(j)=movementData(reachStartIndex(j),timeMarker);
                               tempTitle2 = (strcat( 'starttime:', num2str(t(j))));
                               title(tempTitle2);
                               plot(velocityProfile);
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

                                              fprintf(dataFilepointer2, '%g\t %g \t  %g \t  %g \t %g \n', allFiles,i,j,newStartIndex, newEndIndex);
                                                 
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
                   

                %trialsample(j) = length(trialStartIndex(j):trialEndIndex(j));
                %trialtime (j) = movementData(trialEndIndex(j),timeMarker)-movementData(trialStartIndex(j),timeMarker);
                j = j+1;
                end


                % Calculate Response Time (how long did it take to initiate a response) and Movement Time
                % (how long did it take to reach the target) for each trial             
                 movementStartTime = movementData(reachStartIndex, 3);
                 movementEndTime = movementData(reachEndIndex, 3);
                 RT(1+((i-1) * blockLength):blockLength+((i-1) * blockLength)) = movementStartTime - trialStartTime;
                 MT(1+((i-1) * blockLength):blockLength+((i-1) * blockLength))  = movementEndTime - movementStartTime;            

            end


            cd ..
            cd ..
            cd ..
            close all
            if saveFiles
                backupFileName = [matFileName, num2str(subSeq(numsub)), '.mat'];
                filePath = strcat(pwd,filesep,matFilePath,filesep,'sub_',num2str(subSeq(numsub)),filesep);
                backupFile = [filePath, backupFileName];
                save(backupFile);  
            end

        end
        %clear all;
    end


    ListenChar(0);
    %clear all
    close all
end;               

%Makes a double staircase for bar orientation discrimination (one starts
%at 0 degree difference, the other 10 degree) and after 10 reversals,
%should show us the threshold of accurate orientation discrimination.
%3 down 1 up double staircase, estimate accuracy .792, d' = 1.634
% written by Jianfei, Fall 2015 / modified by Bethany, Summer 2018

% SUBJECTS
% 1 Sydney
% 2 James
% 3 Christian
% 4 Bethany
% 5 Jianfei
% added baseline
% 6 Andrea
% 7 Ryan
% 8 Sydney
% 9 Bethany
% made noise mask
% 10 Jianfei
% 11 Andrea
% 12 Bethany
% 13 Amos

try
    clear all
    VIEWING_DISTANCE_CM = 52;
    MONITOR_WIDTH_CM = 44;
    
    %% subject information %%%%%%%%%%%%%%%%%%%
    KbName('UnifyKeyNames')
    
    subjNum = input('\n Enter subject number: ');
    blockNum = input('\n Enter block number: ');

    oripath = pwd;
    addpath(genpath(oripath));
    %Create a directory for this subject's data if not a practice trial
    if blockNum
        pathdata = strcat(pwd,filesep,'Subject_folders',filesep,'S0',num2str(subjNum),filesep);
        if ~exist(pathdata)
            mkdir(pathdata);
        else
            fprintf('This subject already exists.\n');
        end
        cd(pathdata);
        filenameTxt = strcat(pathdata,filesep,sprintf('%dblock%d',subjNum,blockNum),'_threshold.txt');
        filenameMat = strcat(pathdata,filesep,sprintf('%dblock%d',subjNum,blockNum),'_threshold.mat');
        filenameMatAll = strcat(pathdata,filesep,sprintf('%dblock%d',subjNum,blockNum),'_threshold_all.mat');
    end
    
    %% open window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WhichScreen = max(Screen('Screens'));
    white = WhiteIndex(WhichScreen);
    grey = GrayIndex(WhichScreen);
    
    PsychDebugWindowConfiguration(1,0.5);
    
    Screen('Preference', 'SkipSyncTests', 1);
    [w, winRect] = Screen('OpenWindow',WhichScreen,128);
    if filesep == '\'
        MyCLUT = load('C:\Documents and Settings\js21\My Documents\MATLAB\Bethany\gammaTable1.mat');
        Screen('LoadNormalizedGammaTable', w, MyCLUT.gammaTable1*[1 1 1]);
    end
    [xCen, yCen] = RectCenter(winRect);
    [swidth, ~] = Screen('WindowSize', WhichScreen);
    screenCenter = [xCen yCen];
    screenInfo.center = screenCenter;
    degPerPixel = atan((MONITOR_WIDTH_CM/2)/VIEWING_DISTANCE_CM) * (180/pi) * (2/swidth);
    ppd = 1 / degPerPixel;
    screenInfo.PPD = ppd;
    pixelPerCM = swidth / MONITOR_WIDTH_CM;
    CMperDeg = ppd/pixelPerCM;
    HideCursor;
    secPerFrame = Screen('GetFlipInterval',w);
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    %% staircase settings
    
    nStairDir = 2;
    nHemi = 1;
    nFlankerTilt = 1; % 1 staircase for both
    nCrowdLvl = 2;
    nStaircase = nStairDir*nHemi*nFlankerTilt*nCrowdLvl;
    
    maxori = 20;
    minori = 0.5; % vertical
    ori = repmat([minori maxori],1,nStaircase/2);    % 2 starting points
    oriCatch = repmat([minori maxori],1,nStaircase/2);
    stairdir = repmat([1 0],1,nStaircase/2);  % staircase 1 starts up (1) & 2 starts down (0) direction
    stairdirCatch = repmat([1 0],1,nStaircase/2);  % staircase 1 starts up (1) & 2 starts down (0) direction
    crowdLvl = [0 0 1 1]; % if WhichStair = 1 | 2, no crowding; if WhichStair = 3 | 4, yes crowding
    
    if blockNum % if not practice
        trialNumber = 100*nStaircase; % 100 each for up/down(2)*#nHemi(2)*#flankertilt(2)
    else
        trialNumber = 20;
    end
    
    trials = 0; % initiate trial counter
    flag = 0; % initiate flag 
    
    nReverse = zeros(1,nStaircase);  % counts the number of reversals each staircase
    stairCorrect = zeros(1,nStaircase); % counts # correct in a row each staircase
    trial = zeros(1,nStaircase); % zero trial counters 2 staircases
    realTrial = zeros(1,nStaircase); 
    stimulusReversal = zeros(nStaircase,20); % matrix with stimulus setting at staircase reversals
    
    nReverseCatch = zeros(1,nStaircase);  % counts the number of reversals each staircase
    stairCorrectCatch = zeros(1,nStaircase); % counts # correct in a row each staircase
    stimulusReversalCatch = zeros(nStaircase,20);

    percChange = 0.2; % change signal by 20% at reversals
    
    %% initial value & stimulus settings
    
    key1 = 90; % z = CCW
    key2 = 88; % x = CW
    
    percCatch = 1/3; % 33% catch trials with display in left hemifield
    
    barContrast = .7;
    stimColor = grey - barContrast*grey;
    fixRad = 0.2;  %radius of fixation spot
    barWid = 0.1; 
    barLen = 1.75;
    maskBarWid = 0.15;
    ecc = 14;
    tfDist = 3;
    
    fixRadPx = round(fixRad*ppd);
    FIXATION_POSITION = [xCen - fixRadPx, yCen - fixRadPx, xCen + fixRadPx, yCen + fixRadPx];
    eccPx = round(ecc*ppd);
    tfDistPx = round(tfDist*ppd);
    barLenPx = round(barLen*ppd);
    barWidPx = round(barWid*ppd);
    
    Tmarker = 1;   
    T1 = .1; 
    
    nTotalTrials = ceil(percCatch*trialNumber) + trialNumber;
    hemiIndex(1:trialNumber) = ones; 
    hemiIndex(trialNumber+1:nTotalTrials) = -1;
    x = randperm(nTotalTrials); 
    hemiIndex = hemiIndex(x); 
    idxCatchTrials = find(hemiIndex == -1);
    
    orientList = [-1 1]; % left vs right tilt
    orientIndex = repmat(orientList,1,ceil(nTotalTrials./2));
    n = randperm(nTotalTrials);
    oriIndex = orientIndex(1,n);
    
    flankerTilt = [-1 1];
    flankerIdx = repmat(flankerTilt,1,ceil(nTotalTrials./2));
    n = randperm(nTotalTrials);
    flankerIndex = flankerIdx(1,n);
    
    baseOri = [0 90];
    baseOriIdx = repmat(baseOri,1,ceil(nTotalTrials./2));
    n = randperm(nTotalTrials);
    baseOriIndex = baseOriIdx(1,n);
    
    barVert = ones(barLenPx, barWidPx)*grey - barContrast*grey; % assign pixel values for bar
    barTexVert = Screen('MakeTexture', w, barVert); % convert into texture
    texrect = Screen('Rect', barTexVert); % convert into rect
    
    % destination rects for example bars on instruction screen
    dstRectsInst(:,1) = CenterRectOnPoint(texrect, xCen - eccPx/4, yCen + eccPx/1.5);
    dstRectsInst(:,2) = CenterRectOnPoint(texrect, xCen + eccPx/4, yCen + eccPx/1.5);
    rotAnglesInst = [-20 20];
    
    % initialize stuff for feedback
    stimulus_onset_time(1:nTotalTrials) = zeros;
    r1(1:nTotalTrials) = zeros;
    acc(1:nTotalTrials,1:2) = zeros;
    rt(1:nTotalTrials,1:2) = zeros;
    Keyresponse(1:nTotalTrials,1:2) = zeros;
    stairOrder(1:nTotalTrials) = zeros;
    crowd(1:nTotalTrials) = zeros; 
    rspKey = zeros(1,nTotalTrials);
    rspRatio = [0 0]; % rspRatio(1) counts # lefts, rspRatio(2) counts # rights
    
    %% instructions
    
    KbName('UnifyKeyNames');
    Screen('FillRect', w, grey);
    Screen('TextSize', w, []);
    
    imHolder = [0 0 1000 500];
    [im, ~, ~] = CenterRect(imHolder, winRect); 
    [img, ~, alpha] = imread('screens2.png');
    img(:,:,4) = alpha; 
    imageFinal = Screen('MakeTexture', w, img);
    
    %Practice Instructions
    if blockNum == 0 
        instText = ['Fixate on the center dot during the entire trial.\n'...
            'Press 1 if the center bar is tilted counterclockwise.\n'...
            'Press 2 if the center bar is tilted clockwise.\n\n'...
            'If you respond correctly, you will hear a beep. \n'...
            'If you respond incorrectly, you will hear a lower beep. \n\n'...  
            'If you need to exit the trial, press q. \n\n'...
            'We will start with some practice trials. \n\n'...
            'Press ENTER to continue.'];
    elseif blockNum == 1
        instText = ['The practice block is now complete.\n\n'...
                    'There will be a number of experimental blocks with breaks in between.\n\n'...   
                    'Press ENTER to continue.\n\n'];
    else
        instText = 'Press ENTER to continue.\n\n';
    end
    DrawFormattedText(w, instText, 'Center', 'Center', [255 255 255]);
    
    if blockNum == 0
        Screen('DrawTextures', w, barTexVert, [], dstRectsInst, rotAnglesInst);
    end
    
    Screen('Flip', w);
    FlushEvents('keyDown');
    GetChar;
    Screen('TextSize', w, 20);
    Screen('Flip', w);

    activeKeys = [KbName('1') KbName('2') KbName('q')];
    RestrictKeysForKbCheck(activeKeys);
    
    %% main experiment
    
    while flag == 0 && trials < nTotalTrials % flag = 1 means that we've hit the upper limit
        priorityLevel = MaxPriority(w); % grab high priority to make generating movie as fast as possible
        Priority(priorityLevel);
        
        clear dstRects
        trials = trials + 1 
        WhichStair = randi(nStaircase); % which staircase to use
        stairOrder(trials) = WhichStair; 
        stdir = stairdir(WhichStair);
        crowd(trials) = crowdLvl(WhichStair);
        trial(WhichStair) = trial(WhichStair) + 1;  % count trials on this staircase
        
        if hemiIndex(trials) > 0
            realTrial(WhichStair) = realTrial(WhichStair) + 1;
            stori = ori(WhichStair);
        else
            stori = oriCatch(WhichStair); 
        end
        
        if baseOriIndex(trials)
            cue = {'horizontal','H'};
        else
            cue = {'vertical','V'};
        end
        
        % show fixation cross & pre-cue
        Screen('FillRect', w, grey);
        DrawFormattedText(w, cue{1}, 'Center', FIXATION_POSITION(2)-20, stimColor);
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        Screen('Flip',w);
        WaitSecs(Tmarker);
        
        % draw stimulus display
        dstRects(:,1) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials), yCen); % position the bars; target
        dstRects(:,2) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2)); % UR
        dstRects(:,3) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2)); % LL
        dstRects(:,4) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2)); % UL
        dstRects(:,5) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2)); % LR
        
        if ~crowd(trials)
            dstRects = dstRects(:,1);
        end

        % pick orientations
        r1(trials) = baseOriIndex(trials) + oriIndex(trials)*ori(WhichStair);
        if crowd(trials)
            rotAngles = [r1(trials) repmat(45*flankerIndex(trials),1,4)]; % optimal for threshold elevation
        else
            rotAngles = r1(trials);
        end
        
        DrawFormattedText(w, cue{2}, 'Center', FIXATION_POSITION(2)-20, stimColor);
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles);
        Screen('FillOval', w, stimColor, FIXATION_POSITION, 10);
        Screen('Flip',w);
        stimulus_onset_time(trials) = tic; % Mark the time when the display went up
        WaitSecs(T1);
        
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        Screen('Flip',w);
        
        keypressed = 0;
        while ~keypressed
            [keyIsDown, secs, responseKey] = KbCheck;
            if keyIsDown
                if responseKey(KbName('q'))
                    Screen('CloseAll');
                    ShowCursor;
                end
                if find(responseKey)==key1
                    rspRatio(1) = rspRatio(1) + 1;
                    rspKey(trials) = 0;
                elseif find(responseKey)==key2
                    rspRatio(2) = rspRatio(2) + 1;
                    rspKey(trials) = 1;
                end
                if oriIndex(trials) == 1 && (find(responseKey)==key2) || oriIndex(trials) == -1 && ...
                        (find(responseKey)==key1)
                    acc(trial(WhichStair),WhichStair) = 1;
                else
                    acc(trial(WhichStair),WhichStair) = 0;
                end
                rt(trial(WhichStair),WhichStair) = (secs-stimulus_onset_time(trials));
                Keyresponse(trial(WhichStair),WhichStair) = find(responseKey);
                KbReleaseWait;
                keypressed = 1;
            end
        end
        
        %staircase stuff
        if acc(trial(WhichStair),WhichStair) % IF CORRECT
            Beeper(1000); %play high beep for correct answer
            if hemiIndex(trials) > 0
                stairCorrect(WhichStair) = stairCorrect(WhichStair) + 1;
                if stairCorrect(WhichStair) == 3   % time to make stimulus harder?
                    stairCorrect(WhichStair) = 0;  % zero counter

                    % STAIRCASE DIRECTION FLAGS: 0 = down (getting harder), 1 = up (getting easier)
                    if stairdir(WhichStair) == 1 % if stair direction is currently UP (= 1) this is a reversal
                        stairdir(WhichStair) = 0; % stair direction changes to DOWN
                        nReverse(WhichStair) = nReverse(WhichStair) + 1; % count reversal
                        stimulusReversal(WhichStair, nReverse(WhichStair)) = ori(WhichStair); % record stimulus value ...
                            % (stairStep) at reversal
                    end
                    %change stimulus to make stimulus harder
                    ori(WhichStair) = ori(WhichStair) - percChange*ori(WhichStair); 
                    if ori(WhichStair) < minori
                        ori(WhichStair) = minori;
                    end % can't be negative
                end
            else
                if stairCorrectCatch(WhichStair) == 3   % time to make stimulus harder?
                    stairCorrectCatch(WhichStair) = 0;  % zero counter
                    % STAIRCASE DIRECTION FLAGS: 0 = down (getting harder), 1 = up (getting easier)
                    if stairdirCatch(WhichStair) == 1 % if stair direction is currently UP (= 1) this is a reversal
                        stairdirCatch(WhichStair) = 0; % stair direction changes to DOWN
                        nReverseCatch(WhichStair) = nReverseCatch(WhichStair) + 1; % count reversal
                        stimulusReversalCatch(WhichStair, nReverseCatch(WhichStair)) = oriCatch(WhichStair); % record stimulus value ...
                            % (stairStep) at reversal
                    end
                    %change stimulus to make stimulus harder
                    oriCatch(WhichStair) = oriCatch(WhichStair) - percChange*oriCatch(WhichStair); 
                    if oriCatch(WhichStair) < minori
                        oriCatch(WhichStair) = minori;
                    end % can't be negative
                end
            end
            
        else %incorrect response
            Beeper(500);
            if hemiIndex(trials) > 0
                stairCorrect(WhichStair) = 0;
                if stairdir(WhichStair) == 0
                    stairdir(WhichStair) = 1;   % change direction to UP
                    nReverse(WhichStair) = nReverse(WhichStair) + 1; % count reversal
                    stimulusReversal(WhichStair, nReverse(WhichStair)) = ori(WhichStair); % record stimulus @ this reversal
                end
                %change stimulus to make stimulus easier
                ori(WhichStair) = ori(WhichStair) + percChange*ori(WhichStair);
                if ori(WhichStair) > maxori
                    ori(WhichStair) = maxori;
                end % max of 45
            else
                stairCorrectCatch(WhichStair) = 0;
                if stairdirCatch(WhichStair) == 0
                    stairdirCatch(WhichStair) = 1;   % change direction to UP
                    nReverseCatch(WhichStair) = stairdirCatch(WhichStair) + 1; % count reversal
                    stimulusReversalCatch(WhichStair, stairdirCatch(WhichStair)) = oriCatch(WhichStair); % record stimulus @ this reversal
                end
                %change stimulus to make stimulus easier
                oriCatch(WhichStair) = oriCatch(WhichStair) + percChange*oriCatch(WhichStair);
                if oriCatch(WhichStair) > maxori
                    oriCatch(WhichStair) = maxori;
                end % max of 45
            end
        end
        
        ListenChar(0); %disables keyboard and flushes.
        ListenChar(2); %enables keyboard, no output to command window
        
        % test whether to end experiment
        if (sum(nReverse > 8) == nStaircase) || (sum(realTrial > 100) > 0)
            flag = 1;
        else
            if trials == 200 || trials == 400 || trials == 600
                breakText = 'Please take a break! Press RIGHT ARROW key when ready to resume.\n';
                DrawFormattedText(w, breakText, 'Center', 'Center', [255 255 255]);
                Screen('Flip',w);
                WaitSecs(2); 
                KbWait;
                while KbCheck; end
            end
        end
    end
    
    Screen ('CloseAll');
    ShowCursor;
    
    %--------- DONE --------------------
    % rudimentary data analysis
    if blockNum
        % sum stimulus values at each reversal for separate staircases
%         sumReversal(1) = sum(stimulusReversal(1,4:nReverse(1)));
%         sumReversal(2) = sum(stimulusReversal(2,4:nReverse(2)));
%         
%         % means
%         stair1mean = sumReversal(1)/(nReverse(1)-3);
%         stair2mean = sumReversal(2)/(nReverse(2)-3);
%         
%         % standard deviation
%         StandardDev1 = std(stimulusReversal(1,4:nReverse(1)));
%         StandardDev2 = std(stimulusReversal(2,4:nReverse(2)));
%         
        for i = 1:nStaircase % num staircases 
            final6Rev(:,i) = stimulusReversal(i,nReverse(i)-5:nReverse(i));
        end
        avg1 = mean(final6Rev,1);
        thresholds = [mean(avg1(1:2)) mean(avg1(3:4))];
        
        ListenChar(1); % Turn keyboard output to command window on
        
        dlmwrite(filenameTxt,[WhichStair,trial(WhichStair),stori,r1(trials),acc(trial(WhichStair),WhichStair),...
            stdir,stairCorrect(WhichStair),nReverse(WhichStair),hemiIndex(trials),rspKey(trials),...
            flankerIndex(trials)],'-append','roffset', [],'delimiter', '\t');
        save(filenameTxt);
        save(filenameMatAll);
        save(filenameMat,'trials','trial','r1','acc','nReverse','stimulusReversal','rspRatio','hemiIndex','rspKey',...
            'flankerIndex','stairOrder','realTrial','idxCatchTrial','thresholds');
        
        for i = 1:nStaircase
            plot(stimulusReversal(i,1:nReverse(i)));hold on;
        end
    end
    cd(oripath);
   
catch psychlasterror
    Screen ('CloseAll');
    close all;
    ShowCursor;
    disp(psychlasterror.message);
    disp(psychlasterror.stack);
    
end   % try catch
% end   % function
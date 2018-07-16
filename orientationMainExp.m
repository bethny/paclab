%Makes a double staircase for bar orientation discrimination (one starts
%at 0 degree difference, the other 10 degree) and after 10 reversals,
%should show us the threshold of accurate orientation discrimination.
%3 down 1 up double staircase, estimate accuracy .792, d' = 1.634
% written by Jianfei, Fall 2015 / modified by Bethany, Summer 2018

% CONDITIONS
% 1	same / vert / - flankers
% 2 diff / vert / - flankers
% 3 same / horiz / - flankers
% 4 diff / horiz / - flankers
% 5 same / vert / 0 flankers
% 6 diff / vert / 0 flankers
% 7 same / horiz / 0 flankers
% 8 diff / horiz / 0 flankers
% 9 same / vert / + flankers
% 10 diff / vert / + flankers
% 11 same / horiz / + flankers
% 12 diff / horiz / + flankers

% SUBJECTS
% 1 Bethany
% 2 Jianfei
% 3 Dan
% 4 Bethany
% 5 Sydney
% 6 
% 7 
% 8 
% 9

try
    VIEWING_DISTANCE_CM = 52;
    MONITOR_WIDTH_CM = 44;
    
    %% subject information %%%%%%%%%%%%%%%%%%%
    KbName('UnifyKeyNames')
    
    subjNum = input('\n Enter subject number: ');
    blockNum = input('\n Enter block number: ');
    blockType = input('\n No-action (0) or grasp (1): ');
    task = {'noAction','grasp'};

    oripath = pwd;
%     oripath = 'C:\Documents and Settings\js21\My Documents\MATLAB\Bethany';
    addpath(genpath(oripath));
    cd(oripath);
    %Create a directory for this subject's data if not a practice trial
    pathdata = strcat(pwd,filesep,'Subject_folders',filesep,'S0',num2str(subjNum),filesep);
    if blockNum
        
%         if ~exist(pathdata)
%             mkdir(pathdata);
%         else
%             fprintf('This subject already exists.\n');
%         end
        cd(pathdata);
        filenameTxt = strcat(pathdata,filesep,sprintf('%dblock%d_%s',subjNum,blockNum,task{blockType+1}),'_main.txt');
        filenameMat = strcat(pathdata,filesep,sprintf('%dblock%d_%s',subjNum,blockNum,task{blockType+1}),'_main.mat');
        filenameMatAll = strcat(pathdata,filesep,sprintf('%dblock%d_%s',subjNum,blockNum,task{blockType+1}),'_main_all.mat');
    end
    
    % get threshold info
    data = load(sprintf('%s%dblock1_threshold_all.mat',pathdata,subjNum));
    thresholds = data.thresholds;
    
    %% open window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WhichScreen = max(Screen('Screens'));
    white = WhiteIndex(WhichScreen);
    grey = GrayIndex(WhichScreen);
    
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
    
    %% experiment settings
    
    nTrialPerCnd = 24/2; 
    nFlanker = 3; % +/-/0
    nBaseOri = 2; % 0/90
    nTarg = 2; % same/diff
    nCnd = nFlanker*nTarg*nBaseOri;
    
    trialNumber = nTrialPerCnd*nCnd;
    
    percCatch = 1/3; % 33% catch trials with display in left hemifield
    nCatch = trialNumber*percCatch;
    nTotalTrials = nCatch + trialNumber;
    
    mixedTrialIdx = randperm(nTotalTrials);
    catchIdx = mixedTrialIdx(1:nCatch);
    nonCatchIdx = mixedTrialIdx(nCatch+1:end);
    
    catchCnd = repmat([1:nCnd],1,length(catchIdx)/nCnd);
    x = randperm(nCatch);
    catchCnd = catchCnd(x);
    
    nonCatchCnd = repmat([1:nCnd],1,length(nonCatchIdx)/nCnd);
    x = randperm(nTotalTrials - nCatch);
    nonCatchCnd = nonCatchCnd(x);
    
    cndList = zeros(1,nTotalTrials);
    cndList(catchIdx) = catchCnd;
    cndList(nonCatchIdx) = nonCatchCnd;
    
    targList = [0 1]; % 0 = same, 1 = diff
    baseOriList = [0 90];
    flankerList = [-1 0 1];
    cnds = CombVec(targList,baseOriList,flankerList);
    
    %% initial value & stimulus settings
    
    barContrast = .7;
    stimColor = 1*grey - barContrast*grey;
 
    fixRad = 0.2;  %radius of fixation spot
    barWid = 0.1; 
    barLen = 1.75;
    maskBarWid = 0.15;
    ecc = 14;
    tfDist = 3;
    markerWaitList = [0.75, 1, 1.25];
    mn = 3; % the number of markerWait;
    
    threshMultiplier = 1.5;
    
    T1 = .1; 
    ISI = .3;
    T2 = .1;
    postStim = .05;
  
    hemiIndex = ones(1, nTotalTrials);
    hemiIndex(catchIdx) = -1;
    
    nDiffTrials = sum(mod(cndList,2));
    targTiltList = [-1 1];
    targTiltIdx = repmat(targTiltList,1,nDiffTrials/2);
    n = randperm(nDiffTrials);
    targTiltIdx = targTiltIdx(n);
    
    fixRadPx = round(fixRad*ppd);
    FIXATION_POSITION = [xCen - fixRadPx, yCen - fixRadPx, xCen + fixRadPx, yCen + fixRadPx];
    eccPx = round(ecc*ppd);
    tfDistPx = round(tfDist*ppd);
    barLenPx = round(barLen*ppd);
    barWidPx = round(barWid*ppd);
    
    markerWaitIndex = repmat(markerWaitList,1,floor((nTotalTrials+mn)./mn));
    n = randperm(nTotalTrials);
    markerWait = markerWaitIndex(1,n);
    
    barVert = ones(barLenPx, barWidPx)*grey - barContrast*grey; % assign pixel values for bar
    barTexVert = Screen('MakeTexture', w, barVert); % convert into texture
    texrect = Screen('Rect', barTexVert); % convert into rect
    
    graspSquare = ones(barLenPx, barLenPx)*grey - barContrast*grey;
    graspSquareTex = Screen('MakeTexture', w, graspSquare); % convert into texture
    graspSquareRect = Screen('Rect', graspSquareTex); % convert into rect
    
    % destination rects for example bars on instruction screen
    dstRectsInst(:,1) = CenterRectOnPoint(texrect, xCen - eccPx/4, yCen + eccPx/1.5);
    dstRectsInst(:,2) = CenterRectOnPoint(texrect, xCen + eccPx/4, yCen + eccPx/1.5);
    rotAnglesInst = [-20 20];
    
    % initialize stuff for feedback
    trials = 0; % initiate trial counter
    difftrials = 0;
    stimulus_onset_time(1:nTotalTrials) = zeros;
    targAngle(1:nTotalTrials) = zeros;
    acc(1:nTotalTrials,1:2) = zeros;
    rt(1:nTotalTrials,1:2) = zeros;
    Keyresponse(1:nTotalTrials,1:2) = zeros;
    rspKey = zeros(1,nTotalTrials);
    rspRatio = [0 0]; % rspRatio(1) counts # lefts, rspRatio(2) counts # rights
    
    %% instructions 
    KbName('UnifyKeyNames');
    Screen('FillRect', w, grey);
    Screen('TextSize', w, []);
    
    %Practice Instructions
    if blockNum == 0 
        instText = ['Fixate on the center dot during the entire trial.\n'...
            'Press 1 if the two stimuli are the same.\n'...
            'Press 2 if they are different.\n'...
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
    if ~blockNum
        nTotalTrials = 18;
    end
    while trials < nTotalTrials % flag = 1 means that we've hit the upper limit
        clear dstRects
        clear cDstRects
        
        trials = trials + 1 
        curCnd = cndList(trials); 
        targID(trials) = cnds(1,curCnd);
        baseOri(trials) = cnds(2,curCnd);
        flanker(trials) = cnds(3,curCnd);
        
        if flanker
            threshold = thresholds(2)*threshMultiplier;
        else
            threshold = thresholds(1)*threshMultiplier;
        end
        
        if targID(trials)
            difftrials = difftrials + 1;
            targTilt(trials) = targTiltIdx(difftrials);
        else
            targTilt(trials) = 0;
        end
        
        priorityLevel = MaxPriority(w); % grab high priority to make generating movie as fast as possible
        Priority(priorityLevel);
        
        % show fixation cross
        Screen('FillRect', w, grey);
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        Screen('Flip',w);
        WaitSecs(markerWait(trials));
        
        % draw stimulus display
        dstRects(:,1) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials), yCen); % position the bars; target
        dstRects(:,2) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2));
        dstRects(:,3) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2));
        dstRects(:,4) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2)); % UL
        dstRects(:,5) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2)); % LR
        
        if ~flanker(trials)
            dstRects = dstRects(:,1);
        end

        % pick orientations
        targAngle(trials) = targID(trials)*threshold*targTilt(trials) + baseOri(trials);
        
        if flanker(trials)
            rotAngles1 = [targAngle(trials) repmat(45*flanker(trials),1,4)]; % optimal for threshold elevation
            rotAngles2 = [baseOri(trials) repmat(45*flanker(trials),1,4)];
        else
            rotAngles1 = targAngle(trials);
            rotAngles2 = baseOri(trials);
        end
        
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles1);
        Screen('FillOval', w, stimColor, FIXATION_POSITION, 10);
        Screen('Flip',w);
        stimulus_onset_time(trials) = tic; % Mark the time when the display went up
        WaitSecs(T1);
        
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        Screen('Flip',w);
        WaitSecs(ISI);
        
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles2);
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        Screen('Flip',w);
        WaitSecs(T2);
        
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        Screen('Flip',w);
        
        tic
        % NOISE CIRCLE
        circRect = SetRect(0,0, barLenPx, barLenPx);
        cDstRects(:,1) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials), yCen); % position the bars; target
        if flanker(trials)
            cDstRects(:,2) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), ...
                yCen - tfDistPx/sqrt(2));
            cDstRects(:,3) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), ...
                yCen + tfDistPx/sqrt(2));
            cDstRects(:,4) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), ...
                yCen - tfDistPx/sqrt(2));
            cDstRects(:,5) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), ...
                yCen + tfDistPx/sqrt(2));
        end
        circAperture = Screen('OpenOffscreenwindow', w, 128, circRect); % offscreen aperture for circle mask
        sqAperture = Screen('OpenOffscreenwindow', w, 128, circRect);
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        AssertOpenGL;
        contrast = 2;
        tex = CreateProceduralNoise(w, barLenPx, barLenPx, 'Perlin', [0.0 0.0 0.0 0.0]); % 0.5 0.5 0.5 0 for grey
        Screen('DrawTextures', w, tex, [], cDstRects, [], 0, [], [1 1 1], [], [], repmat([contrast, 0, 0, 0],5,1)');
        Screen('FillRect', sqAperture, [255 255 255 0], circRect); % alpha = 0 so noise can come thru
        Screen('DrawTextures', w, sqAperture, [], cDstRects(:,1), [], 0); 
        if flanker(trials)
            Screen('FillOval', circAperture, [255 255 255 0], circRect); % alpha = 0 so noise can come thru
            Screen('DrawTextures', w, circAperture, [], cDstRects(:,2:5), [], 0);        
        end
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);

        diff = T2 + postStim - toc;
        WaitSecs(diff);
        Screen('Flip',w);
        
        keypressed = 0;
        while ~keypressed
            [keyIsDown, secs, responseKey] = KbCheck;
            if keyIsDown
                if responseKey(KbName('q'))
                    Screen('CloseAll');
                    ShowCursor;
                end
                if responseKey(KbName('1'))
                    rspKey(trials) = 0;
                elseif responseKey(KbName('2'))
                    rspKey(trials) = 1;
                end
                if (~targID(trials) && responseKey(KbName('1'))) || (targID(trials) && responseKey(KbName('2'))) % 1 for same, 2 for diff
                    acc(trials) = 1;
                    Beeper(1000);
                else
                    acc(trials) = 0;
                    Beeper(500);
                end
                rt(trials) = (secs-stimulus_onset_time(trials));
                Keyresponse(trials) = find(responseKey);
                KbReleaseWait;
                keypressed = 1;
            end
        end
        
        ListenChar(0); %disables keyboard and flushes.
        ListenChar(2); %enables keyboard, no output to command window
        
        % test whether to offer a break
        if trials == 48 || trials == 96 || trials == 144
            breakText = 'Please take a break! Press 1 when ready to resume.\n';
            DrawFormattedText(w, breakText, 'Center', 'Center', [255 255 255]);
            Screen('Flip',w);
            WaitSecs(2); 
            KbWait;
            while KbCheck; end
        end
        end
    
    Screen ('CloseAll');
    ShowCursor;
    
    %--------- DONE --------------------
    % rudimentary data analysis
    if blockNum
        ListenChar(1); % Turn keyboard output to command window on
        
%         dlmwrite(filenameTxt,[WhichStair,trial(WhichStair),stori,targAngle(trials),acc(trial(WhichStair),WhichStair),...
%             stdir,stairCorrect(WhichStair),nReverse(WhichStair),hemiIndex(trials),rspKey(trials),...
%             flankerIndex(trials)],'-append', ...
%             'roffset', [],'delimiter', '\t');
%         save(filenameTxt);
        save(filenameMatAll);
        save(filenameMat,'trials','acc','rspRatio','hemiIndex','rspKey','cndList','targTilt','baseOri',...
            'flanker','targAngle','targID');
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
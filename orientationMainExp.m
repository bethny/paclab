%Makes a double staircase for bar orientation discrimination (one starts
%at 0 degree difference, the other 10 degree) and after 10 reversals,
%should show us the threshold of accurate orientation discrimination.
%3 down 1 up double staircase, estimate accuracy .792, d' = 1.634
% written by Jianfei, Fall 2015 / modified by Bethany, Summer 2018

% CONDITIONS
% 1	ori = -, base = vert, crowding = -
% 2 ori = 0, base = vert, crowding = -
% 3 ori = +, base = vert, crowding = -
% 4 ori = -, base = horiz, crowding = -
% 5 ori = 0, base = horiz, crowding = -
% 6 ori = +, base = horiz, crowding = -
% 7 ori = -, base = vert, crowding = 0
% 8 ori = 0, base = vert, crowding = 0
% 9 ori = +, base = vert, crowding = 0
% 10 ori = -, base = horiz, crowding = 0
% 11 ori = 0, base = horiz, crowding = 0
% 12 ori = +, base = horiz, crowding = 0
% 13 ori = -, base = vert, crowding = +
% 14 ori = 0, base = vert, crowding = +
% 15 ori = +, base = vert, crowding = +
% 16 ori = -, base = horiz, crowding = +
% 17 ori = 0, base = horiz, crowding = +
% 18 ori = +, base = horiz, crowding = +

% SUBJECTS
% 1 
% 2 
% 3 
% 4 
% 5 
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
    
    nTilt = 3; % +/-/0
    nBaseOri = 2; % 0/90
    nCrwd = 3; % +45/-45/none
    nCnd = nTilt*nCrwd*nBaseOri;
    
    if blockNum % if not practice
        trialNumber = 40*nCnd;
    else
        trialNumber = 20;
    end
    
    percCatch = .35; % 35% catch trials with display in left hemifield
    nTotalTrials = ceil(percCatch*trialNumber) + trialNumber;
    
    cndList = repmat([1:nCnd],1,ceil(nTotalTrials/nCnd));
    x = randperm(nTotalTrials);
    cndList = cndList(x);
    
    threshold = 8;
    none = 0; % vertical
    
    tiltList = [-1*threshold 0 threshold];
    baseOriList = [0 90];
    crowdList = [-1 0 1];
    cnds = combvec(tiltList,baseOriList,crowdList);
    
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
    
    T1 = .1; 
    ISI = .15;
    T2 = .1;
    postStim = .05;
    graspDur = .6;
  
    hemiIndex(1:trialNumber) = ones; 
    hemiIndex(trialNumber+1:nTotalTrials) = -1;
    x = randperm(nTotalTrials); 
    hemiIndex = hemiIndex(x); 
    
%     orientList = [-1 1]; % left vs right tilt
%     orientIndex = repmat(orientList,1,ceil(nTotalTrials./2));
%     n = randperm(nTotalTrials);
%     oriIndex = orientIndex(1,n);
%     
%     baseOriList = [0 90]; % vertical vs horizontal bar
%     baseOriIdx = repmat(baseOriList,1,ceil(nTotalTrials./2));
%     n = randperm(nTotalTrials);
%     baseOriIndex = baseOriIdx(1,n);
%     
%     flankerTilt = [-1 1];
%     flankerIdx = repmat(flankerTilt,1,ceil(nTotalTrials./2));
%     n = randperm(nTotalTrials);
%     flankerIndex = flankerIdx(1,n);
    
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
    
    % destination rects for example bars on instruction screen
    dstRectsInst(:,1) = CenterRectOnPoint(texrect, xCen - eccPx/4, yCen + eccPx/1.5);
    dstRectsInst(:,2) = CenterRectOnPoint(texrect, xCen + eccPx/4, yCen + eccPx/1.5);
    rotAnglesInst = [-20 20];
    
    % initialize stuff for feedback
    trials = 0; % initiate trial counter
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
    
    imHolder = [0 0 1000 500];
    [im, ~, ~] = CenterRect(imHolder, winRect); 
    [img, ~, alpha] = imread('screens2.png');
    img(:,:,4) = alpha; 
    imageFinal = Screen('MakeTexture', w, img);
    
    %Practice Instructions
    if blockNum == 0 
        instText = ['Fixate on the center dot during the entire trial.\n'...
            'Press the LEFT ARROW KEY if the center bar\n'...
            'is tilted to the left (counterclockwise).\n'...
            'Press the RIGHT ARROW KEY if the center bar\n'...
            'is tilted to the right (clockwise).\n\n'...
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
    
%     if blockNum == 0
%         Screen('DrawTexture', w, imageFinal, [], im);
%         Screen('Flip', w);
%         FlushEvents('keyDown');
%         GetChar;
%         Screen('TextSize', w, 20);
%         Screen('Flip', w);
%     end

    activeKeys = [KbName('1') KbName('2') KbName('q')];
    RestrictKeysForKbCheck(activeKeys);
    
    %% main experiment
    
    while trials < nTotalTrials % flag = 1 means that we've hit the upper limit
        clear dstRects
        clear cDstRects
        
        trials = trials + 1 
        curCnd = cndList(trials); 
        targOri = cnds(1,curCnd);
        baseOri = cnds(2,curCnd);
        crowding = cnds(3,curCnd);
        
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
        
        if ~crowd(trials)
            dstRects = dstRects(:,1);
        end

        % pick orientations
        targAngle = targOri + baseOri;
        
        if crowding
            rotAngles1 = [targAngle repmat(45*crowding,1,4)]; % optimal for threshold elevation
            rotAngles2 = [baseOri repmat(45*crowding,1,4)];
        else
            rotAngles1 = targAngle;
            rotAngles2 = baseOri;
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
        if crowd(trials)
            cDstRects(:,2) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), ...
                yCen - tfDistPx/sqrt(2));
            cDstRects(:,3) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), ...
                yCen + tfDistPx/sqrt(2));
            cDstRects(:,4) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), ...
                yCen - tfDistPx/sqrt(2));
            cDstRects(:,5) = CenterRectOnPoint(circRect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), ...
                yCen + tfDistPx/sqrt(2));
        end
        aperture = Screen('OpenOffscreenwindow', w, 128, circRect); % offscreen aperture for circle mask
        Screen('FillOval', aperture, [255 255 255 0], circRect); % alpha = 0 so noise can come thru
        Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        AssertOpenGL;
        contrast = 2;
        tex = CreateProceduralNoise(w, barLenPx, barLenPx, 'Perlin', [0.0 0.0 0.0 0.0]); % 0.5 0.5 0.5 0 for grey
        Screen('DrawTextures', w, tex, [], cDstRects, [], 0, [], [1 1 1], [], [], repmat([contrast, 0, 0, 0],3,1)');
        Screen('DrawTextures', w, aperture, [], cDstRects, [], 0);        
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);

        diff = T2 + postStim - toc;
        WaitSecs(diff);
        Screen('Flip',w);
        WaitSecs(graspDur);
        
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
                if responseKey(KbName('1'))
                    rspKey(trials) = 0;
                elseif responseKey(KbName('2'))
                    rspKey(trials) = 1;
                end
                if ~targOri && responseKey(KbName('1')) || targOri && responseKey(KbName('2')) % 1 for same, 2 for diff
                    acc(trials) = 1;
                else
                    acc(trials) = 0;
                end
                rt(trials) = (secs-stimulus_onset_time(trials));
                Keyresponse(trials) = find(responseKey);
                KbReleaseWait;
                keypressed = 1;
            end
        end
        
        %staircase stuff
        if ((find(responseKey)==key1)&&(targAngle==0)) || ((find(responseKey)==key2)&&(targAngle~=0))
            acc(trial(WhichStair),WhichStair) = 1; % IF CORRECT
            Beeper(1000); %play high beep for correct answer
        else %incorrect response
            Beeper(500);
        end
        
        if hemiIndex > 0
            realTrials = realTrials + 1;
        end
        
        ListenChar(0); %disables keyboard and flushes.
        ListenChar(2); %enables keyboard, no output to command window
        
        % test whether to offer a break
        if trials == 200 || trials == 400 || trials == 600
            breakText = 'Please take a break! Press RIGHT ARROW key when ready to resume.\n';
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
        % sum stimulus values at each reversal for separate staircases
        sumReversal(1) = sum(stimulusReversal(1,4:nReverse(1)));
        sumReversal(2) = sum(stimulusReversal(2,4:nReverse(2)));
        
        % means
        stair1mean = sumReversal(1)/(nReverse(1)-3);
        stair2mean = sumReversal(2)/(nReverse(2)-3);
        
        % standard deviation
        StandardDev1 = std(stimulusReversal(1,4:nReverse(1)));
        StandardDev2 = std(stimulusReversal(2,4:nReverse(2)));
        
        ListenChar(1); % Turn keyboard output to command window on
        
        dlmwrite(filenameTxt,[WhichStair,trial(WhichStair),stori,targAngle(trials),acc(trial(WhichStair),WhichStair),...
            stdir,stairCorrect(WhichStair),nReverse(WhichStair),hemiIndex(trials),rspKey(trials),...
            flankerIndex(trials)],'-append', ...
            'roffset', [],'delimiter', '\t');
        save(filenameTxt);
        save(filenameMatAll);
        save(filenameMat,'trials','trial','r1','acc','nReverse','stimulusReversal','rspRatio','hemiIndex','rspKey',...
            'flankerIndex','stairOrder','realTrial');
        
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
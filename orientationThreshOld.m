%Makes a double staircase for bar orientation discrimination (one starts
%at 0 degree difference, the other 10 degree) and after 10 reversals,
%should show us the threshold of accurate orientation discrimination.
%3 down 1 up double staircase, estimate accuracy .792, d' = 1.634
% written by Jianfei, Fall 2015 / modified by Bethany, Summer 2018

% make it thinner, shorter, smaller flanker-target distance
% see if you need a pure vertical line -- literature review 
% maybe grasp placeholder should be the circle with a + inside

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

try
    VIEWING_DISTANCE_CM = 52;
    MONITOR_WIDTH_CM = 44;
    
    %% subject information %%%%%%%%%%%%%%%%%%%
    KbName('UnifyKeyNames')
    subjNum = input('\n Enter subject number: ');
    blockNum = input('\n Enter block number: ');
    nStim = input('\n T1 and T2? 0 for no, 1 for yes: ');
    lowContrast = input('\n Low contrast? 0 for no, 1 for yes: ');
    oripath = pwd;
    addpath(genpath(oripath));
    %Create a directory for this subject's data if not a practice trial
    if blockNum
        pathdata = strcat(pwd,filesep,'Subject_folders',filesep,num2str(subjNum),filesep);
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
    [swidth, sheight]=Screen('WindowSize', WhichScreen);
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
    
    %% staircase value
    
    nStairDir = 2;
    nHemi = 2;
%     nFlankerTilt = 2;
    nFlankerTilt = 1; % 1 staircase for both
    nCrowdLvl = 2;
    nStaircase = nStairDir*nHemi*nFlankerTilt*nCrowdLvl;
    
    if blockNum % if not practice
        trialNumber = 100*nStairDir*nHemi*nFlankerTilt*nCrowdLvl; % 100 each for up/down(2)*#nHemi(2)*#flankertilt(2)
    else
        trialNumber = 20;
    end
    
    trials = 0; % initiate trial counter
    stairdir = repmat([1 0],1,nHemi*nFlankerTilt*nCrowdLvl);  % staircase 1 starts up (1) & 2 starts down (0) direction
    nReverse = zeros(1,nStaircase);  % counts the number of reversals each staircase
    stairCorrect = zeros(1,nStaircase); % counts # correct in a row each staircase
    trial = zeros(1,nStaircase); % zero trial counters 2 staircases
    stimulusReversal = zeros(nStaircase,20); % matrix with stimulus setting at staircase reversals
    rspKey = zeros(1,trialNumber);
    rspRatio = [0 0]; % rspRatio(1) counts # lefts, rspRatio(2) counts # rights
    flag = 0;
    maxori = 20;
    minori = 0.5; % vertical
    ori = repmat([minori maxori],1,nHemi*nFlankerTilt*nCrowdLvl);    % 2 starting points
    delta = 0.5;          % how much to change signal at reversals
    
    presentationDur = .1;
    
    T1 = .1; 
    ISI = .2;
    T2 = .1;
    
    %% initial value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    barContrast = .7;
    
    fixRad = 0.2;  %radius of fixation spot
%     barWid = 0.3; % from bulakowski
    barWid = 0.1; 
%     barLen = 2.5; % from bulakowski
    barLen = 1.75;
    maskBarWid = 0.15;
    ecc = 14;
%     tfDist = 4; % distance between target + flanker
    tfDist = 3;
    markerWaitList = [0.75, 1, 1.25];
    mn = 3;     % the number of markerWait;
    
    orientList = [-1 1]; % left vs right tilt
    orientIndex = repmat(orientList,1,trialNumber./2);
    n = randperm(trialNumber);
    oriIndex = orientIndex(1,n);
    
%     hemiList = [-1 1]; % left vs right hemi
%     hemiIdx = repmat(hemiList,1,trialNumber./2);
%     n = randperm(trialNumber);
%     hemiIndex = hemiIdx(1,n);
    
    flankerTilt = [-1 1];
    flankerIdx = repmat(flankerTilt,1,trialNumber./2);
    n = randperm(trialNumber);
    flankerIndex = flankerIdx(1,n);
    
    fixRadPx = round(fixRad*ppd);
    FIXATION_POSITION = [xCen - fixRadPx, yCen - fixRadPx, xCen + fixRadPx, yCen + fixRadPx];
    eccPx = round(ecc*ppd);
    tfDistPx = round(tfDist*ppd);
    barLenPx = round(barLen*ppd);
    barWidPx = round(barWid*ppd);
    maskBarWidPx = round(maskBarWid*ppd);
    
    markerWaitIndex = repmat(markerWaitList,1,floor((trialNumber+mn)./mn));
    n = randperm(trialNumber);
    markerWait = markerWaitIndex(1,n);
    
    if lowContrast
        barVert = ones(barLenPx, barWidPx)*grey - barContrast*grey; % assign pixel values for bar
    else
        barVert = ones(barLenPx, barWidPx);
    end
    barTexVert = Screen('MakeTexture', w, barVert); % convert into texture
    texrect = Screen('Rect', barTexVert); % convert into rect
    
    circrect = [0 0 barLenPx barLenPx]; 
    
    dstRectsInst(:,1) = CenterRectOnPoint(texrect, xCen - eccPx/4, yCen + eccPx/1.5);
    dstRectsInst(:,2) = CenterRectOnPoint(texrect, xCen + eccPx/4, yCen + eccPx/1.5);
    rotAnglesInst = [-20 20];
    
    %initialize stuff for feedback
    stimulus_onset_time(1:trialNumber) = zeros;
    r1(1:trialNumber) = zeros;
    acc(1:trialNumber,1:2) = zeros;
    rt(1:trialNumber,1:2) = zeros;
    Keyresponse(1:trialNumber,1:2) = zeros;
    hemiIndex(1:trialNumber) = zeros;
    if nFlankerTilt == 2
        flankerIndex(1:trialNumber) = zeros;
    end
    stairOrder(1:trialNumber) = zeros;
    crowd(1:trialNumber) = zeros; 
    
    %% instructions
    
    KbName('UnifyKeyNames');
    Screen('FillRect', w, grey);
    Screen('TextSize', w, []);
    
%     imHolder = [0 0 1000 500];
%     [im, ~, ~] = CenterRect(imHolder, winRect); 
%     [img, ~, alpha] = imread('screens2.png');
%     img(:,:,4) = alpha; 
%     imageFinal = Screen('MakeTexture', w, img);
    
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
%     
%     if blockNum == 0 && ~nStim
%         Screen('DrawTexture', w, imageFinal, [], im);
%         Screen('Flip', w);
%         FlushEvents('keyDown');
%         GetChar;
%         Screen('TextSize', w, 20);
%         Screen('Flip', w);
%     end

    activeKeys = [KbName('LeftArrow') KbName('RightArrow') KbName('q')];
    RestrictKeysForKbCheck(activeKeys);
    
    %% main experiment
    
    while flag == 0 && trials < trialNumber % flag = 1 means that we've hit the upper limit
        clear dstRects
        trials = trials + 1 
        WhichStair = randi(nStaircase); % which staircase to use
        stairOrder(trials) = WhichStair; 
        stdir = stairdir(WhichStair);
        stori = ori(WhichStair);
        trial(WhichStair) = trial(WhichStair) + 1;  % count trials on this staircase
        
        if WhichStair == 1 || WhichStair == 2
            hemiIndex(trials) = -1;
            if nFlankerTilt == 2
                flankerIndex(trials) = -1;
            end
            if nCrowdLvl == 2
                crowd = 0;
            end
        elseif WhichStair == 3 || WhichStair == 4
            hemiIndex(trials) = -1;
            if nFlankerTilt == 2
                flankerIndex(trials) = 1;
            end
            if nCrowdLvl == 2
                crowd = 1;
            end
        elseif WhichStair == 5 || WhichStair == 6
            hemiIndex(trials) = 1;
            if nFlankerTilt == 2
                flankerIndex(trials) = -1;
            end
            if nCrowdLvl == 2
                crowd = 0;
            end
        else
            hemiIndex(trials) = 1;
            if nFlankerTilt == 2
                flankerIndex(trials) = 1;
            end
            if nCrowdLvl == 2
                crowd = 1;
            end
        end
        
        priorityLevel = MaxPriority(w); % grab high priority to make generating movie as fast as possible
        Priority(priorityLevel);
        
        % show fixation cross
        Screen('FillRect', w, grey);
        % Draw fixation to indicate the start of the trial
        Screen('FillOval', w, 1*grey - barContrast*grey,FIXATION_POSITION,10);
        Screen('Flip',w);
        WaitSecs(markerWait(trials));
        
        % draw stimulus display
        % position the stimuli 
        dstRects(:,1) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials), yCen); % position the bars; target
        dstRects(:,2) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2));
        dstRects(:,3) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2));
        
        upperBound = dstRects(2,2);
        lowerBound = dstRects(4,3);
        leftBound = dstRects(1,3);
        rightBound = dstRects(3,2);
        
        if ~crowd
            dstRects = dstRects(:,1);
        end
            
        % generate & position the mask
        nRow = 4;
        nCol = 4;  
        point = gridGen(nRow, nCol, upperBound, lowerBound, leftBound, rightBound);
        nStimBar = 3;
        maskBarLen = (barLen*barWid*nStimBar)/(maskBarWid*length(point)); % you want the same nun of pixels on screen
        maskBarLenPx = round(maskBarLen*ppd);
        if lowContrast
            maskBar = ones(maskBarLenPx, maskBarWidPx)*grey - barContrast*grey;
        else
            maskBar = ones(maskBarLenPx, maskBarWidPx);
        end
        maskBarTex = Screen('MakeTexture', w, maskBar);
        maskTexrect = Screen('Rect', maskBarTex); % convert into rect
        
        % position the circle mask
        circDstRects = CenterRectOnPoint(circrect, xCen + eccPx*hemiIndex(trials), yCen); 
        circColor = [0 0 0 150]; 
        
        % pick orientations
        r1(trials) = oriIndex(trials)*ori(WhichStair);
        if crowd
            rotAngles = [r1(trials) 45*flankerIndex(trials) 45*flankerIndex(trials)]; % optimal for threshold elevation
            rotAnglesT2 = [0 rotAngles(2:3)];
        else
            rotAngles = r1(trials);
            rotAnglesT2 = 0;
        end
        
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles);
%         Screen('DrawText', w, sprintf('%g, stair = %d, correct = %d',r1(trials),WhichStair,stairCorrect(WhichStair)), ...
%             30, 30, [255, 255, 255]);
        Screen('FillOval', w, 1*grey - barContrast*grey, FIXATION_POSITION, 10);
        Screen('Flip',w);
        stimulus_onset_time(trials) = tic; % Mark the time when the display went up
               
        % draw mask display
        if ~nStim
            WaitSecs(presentationDur);
            
            if ~lowContrast
                for i = 1:length(point)
                    maskDstRects(:,i) = CenterRectOnPoint(maskTexrect, point(i,1)+randi([-30 30]), point(i,2)+randi([-30 30]));
                end
                maskRotAngles = randi(360,1,length(point));
                Screen('FillOval', w, white,FIXATION_POSITION,10);
                Screen('DrawTextures', w, maskBarTex, [], maskDstRects, maskRotAngles);
                Screen('Flip',w);
                WaitSecs(maskDur);
            end
        else
            WaitSecs(T1);
            
            Screen('FillOval', w, 1*grey - barContrast*grey,FIXATION_POSITION,10);
            Screen('Flip',w);
            WaitSecs(ISI);
            
            Screen('DrawTextures', w, barTexVert, [], dstRects, rotAnglesT2);
            Screen('FillOval', w, 1*grey - barContrast*grey,FIXATION_POSITION,10);
            Screen('Flip',w);
            WaitSecs(T2);
        end
        
        Screen('FillOval', w, 1*grey - barContrast*grey,FIXATION_POSITION,10);
        Screen('FrameOval', w, 1*grey - barContrast*grey, circDstRects);
        Screen('DrawTextures', w, barTexVert, [], dstRects(:,1), 0);
        Screen('DrawTextures', w, barTexVert, [], dstRects(:,1), 90);
        Screen('Flip',w);
        
        keypressed = 0;
        while ~keypressed
            [keyIsDown, secs, responseKey] = KbCheck;
            if keyIsDown
                if responseKey(KbName('q'))
                    Screen('CloseAll');
                    ShowCursor;
                end
                if responseKey(KbName('LeftArrow'))
                    rspRatio(1) = rspRatio(1) + 1;
                    rspKey(trials) = 0;
                elseif responseKey(KbName('RightArrow'))
                    rspRatio(2) = rspRatio(2) + 1;
                    rspKey(trials) = 1;
                end
                if oriIndex(trials) == 1 && responseKey(KbName('RightArrow')) || oriIndex(trials) == -1 && ...
                        responseKey(KbName('LeftArrow'))
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
        percChange = 0.2; % 20%, used by Livne/Sagi
        if acc(trial(WhichStair),WhichStair) % IF CORRECT
            Beeper(1000); %play high beep for correct answer
            stairCorrect(WhichStair) = stairCorrect(WhichStair) + 1;
            if stairCorrect(WhichStair) == 3   % time to make stimulus harder?
                stairCorrect(WhichStair) = 0;  % zero counter
                
                % STAIRCASE DIRECTION FLAGS: 0 = down (getting harder), 1 = up (getting easier)
                if stairdir(WhichStair) == 1 % if stair direction is currently UP (= 1) this is a reversal
                    stairdir(WhichStair) = 0; % stair direction changes to DOWN
                    nReverse(WhichStair) = nReverse(WhichStair) + 1; % count reversal
                    stimulusReversal(WhichStair, nReverse(WhichStair)) = ori(WhichStair); % record stimulus value (stairStep)                                                                         % at reversal
                end
                %change stimulus to make stimulus harder
                ori(WhichStair) = ori(WhichStair) - percChange*ori(WhichStair); 
                if ori(WhichStair) < minori
                    ori(WhichStair) = minori;
                end % can't be negative
            end
            
        else %incorrect response
            Beeper(500);
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
        end
        
        ListenChar(0); %disables keyboard and flushes.
        ListenChar(2); %enables keyboard, no output to command window
        
        % test whether to end experiment
        if (nReverse(1) > 8 &&  nReverse(2) > 8 &&  nReverse(3) > 8 &&  nReverse(4) > 8 &&  nReverse(5) > 8 ...
             &&  nReverse(6) > 8 &&  nReverse(7) > 8 &&  nReverse(8) > 8) || trial(1) > 100 || trial(2) > 100 ...
             || trial(3) > 100 || trial(4) > 100 || trial(5) > 100 || trial(6) > 100 || trial(7) > 100 ...
             || trial(8) > 100
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
        sumReversal(1) = sum(stimulusReversal(1,4:nReverse(1)));
        sumReversal(2) = sum(stimulusReversal(2,4:nReverse(2)));
        
        % means
        stair1mean = sumReversal(1)/(nReverse(1)-3);
        stair2mean = sumReversal(2)/(nReverse(2)-3);
        
        % standard deviation
        StandardDev1 = std(stimulusReversal(1,4:nReverse(1)));
        StandardDev2 = std(stimulusReversal(2,4:nReverse(2)));
        
        ListenChar(1); % Turn keyboard output to command window on
        
        dlmwrite(filenameTxt,[WhichStair,trial(WhichStair),stori,r1(trials),acc(trial(WhichStair),WhichStair),...
            stdir,stairCorrect(WhichStair),nReverse(WhichStair),hemiIndex(trials),rspKey(trials),...
            flankerIndex(trials)],'-append', ...
            'roffset', [],'delimiter', '\t');
        save(filenameTxt);
        save(filenameMatAll);
        save(filenameMat,'trials','trial','r1','acc','nReverse','stimulusReversal','rspRatio','hemiIndex','rspKey',...
            'flankerIndex','stairOrder');
        
        for i = 1:8
            plot(stimulusReversal(i,1:nReverse(i)));hold on;
        end
    end
    cd(oripath);
    
%     saveas(gcf, pathdata, 'fig');
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 10 10];
%     print(pathdata,'-dpng')
%     fprintf('Figure saved.\n')
%     
catch psychlasterror
    Screen ('CloseAll');
    close all;
    ShowCursor;
    disp(psychlasterror.message);
    disp(psychlasterror.stack);
    %psychrethrow (psychlasterror);
    
end   % try catch
% end   % function
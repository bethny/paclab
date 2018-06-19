%Makes a double staircase for bar orientation discrimination (one starts
%at 0 degree difference, the other 10 degree) and after 10 reversals,
%should show us the threshold of accurate orientation discrimination.
%3 down 1 up double staircase, estimate accuracy .792, d' = 1.634
% written by Jianfei, Fall 2015 / modified by Bethany, Summer 2018

% SUBJECTS
% BLOCKING: A = right side only, B = both hemifields
% 1 AB Jianfei
% 2 BA Bethany
% 3 AB Christian
% 4 BA

try
    VIEWING_DISTANCE_CM = 52;
    MONITOR_WIDTH_CM = 44;
    
    %% subject information %%%%%%%%%%%%%%%%%%%
    KbName('UnifyKeyNames')
    subjNum = input('\n Enter subject number: ');
    blockNum = input('\n Enter block number: ');
    hemifield = input('\n Num hemifields? 1 or 2: ');
    oripath = pwd;
%     addpath(genpath('C:\Documents and Settings\js21\My Documents\MATLAB\Bethany\')); 
    addpath(genpath(oripath));
    %Create a directory for this subject's data if not a practice trial
    if blockNum
        pathdata = strcat(pwd,filesep,'Subject_folders',filesep,num2str(subjNum),filesep);
        if ~exist(pathdata)
            mkdir(pathdata);
        end
        cd(pathdata);
        filenameTxt = strcat(pathdata,filesep,sprintf('%dblock%d_%dhemi',subjNum,blockNum,hemifield),'_threshold.txt');
        filenameMat = strcat(pathdata,filesep,sprintf('%dblock%d_%dhemi',subjNum,blockNum,hemifield),'_threshold.mat');
        filenameMatAll = strcat(pathdata,filesep,sprintf('%dblock%d_%dhemi',subjNum,blockNum,hemifield),'_threshold_all.mat');
    end
    
    %% open window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WhichScreen = max(Screen('Screens'));
    white = WhiteIndex(WhichScreen);
    grey = GrayIndex(WhichScreen);
    
%     Screen('Preference', 'SkipSyncTests', 1);
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
    
    %% staircase value
    if blockNum % if not practice
        trialNumber = 100*hemifield*2;
%         trialNumber = 40;
    else
        trialNumber = 20;
    end
       
    trials = 0; % initiate trial counter
    stairdir = repmat([1 0],1,hemifield);  % staircase 1 starts up (1) & 2 starts down (0) direction
    nReverse = zeros(1,hemifield*2);  % counts the number of reversals each staircase
    stairCorrect = zeros(1,hemifield*2); % counts # correct in a row each staircase
    trial = zeros(1,hemifield*2); % zero trial counters 2 staircases
    stimulusReversal = zeros(hemifield*2,20); % matrix with stimulus setting at staircase reversals
    rspRatio = [0 0]; % rspRatio(1) counts # lefts, rspRatio(2) counts # rights
    flag = 0;
    maxori = 10;
    minori = 0.5; % vertical
    ori = repmat([minori maxori],1,hemifield);    % 2 starting points
    delta = 0.5;          % how much to change signal at reversals\
    presentationDur = .2;
    maskDur = .2;
    
    %% initial value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fixRad = 0.2;  %radius of fixation spot
    barWid = 0.3; % from bulakowski
    barLen = 2.5; % from bulakowski
    maskBarWid = 0.15;
    maskBarLen = 1.5;
    ecc = 14;
    tfDist = 4; % distance between target + flanker
    markerWaitList = [0.75, 1, 1.25];
    mn = 3;     % the number of markerWait;
    
    orientList = [-1 1]; % left vs right tilt
    orientIndex = repmat(orientList,1,trialNumber./2);
    n = randperm(trialNumber);
    oriIndex = orientIndex(1,n);
    
    if hemifield == 2
        hemiList = [-1 1]; % left vs right hemifield
    else
        hemiList = [1 1];
    end
    hemiIdx = repmat(hemiList,1,trialNumber./2);
    n = randperm(trialNumber);
    hemiIndex = hemiIdx(1,n);
    
    fixRadPx = round(fixRad*ppd);
    FIXATION_POSITION = [xCen - fixRadPx, yCen - fixRadPx, xCen + fixRadPx, yCen + fixRadPx];
    eccPx = round(ecc*ppd);
    tfDistPx = round(tfDist*ppd);
    barLenPx = round(barLen*ppd);
    barWidPx = round(barWid*ppd);
    maskBarWidPx = round(maskBarWid*ppd);
    maskBarLenPx = round(maskBarLen*ppd);
    
    markerWaitIndex = repmat(markerWaitList,1,floor((trialNumber+mn)./mn));
    n = randperm(trialNumber);
    markerWait = markerWaitIndex(1,n);
    
    barVert = ones(barLenPx, barWidPx); % assign pixel values for bar
    barTexVert = Screen('MakeTexture', w, barVert); % convert into texture
    texrect = Screen('Rect', barTexVert); % convert into rect
    
    dstRectsInst(:,1) = CenterRectOnPoint(texrect, xCen - eccPx/4, yCen + eccPx/1.5);
    dstRectsInst(:,2) = CenterRectOnPoint(texrect, xCen + eccPx/4, yCen + eccPx/1.5);
    rotAnglesInst = [-20 20];
    
    maskBar = ones(maskBarLenPx, maskBarWidPx);
    maskBarTex = Screen('MakeTexture', w, maskBar);
    maskTexrect = Screen('Rect', maskBarTex); % convert into rect
    
    %initialize stuff for feedback
    stimulus_onset_time(1:trialNumber) = zeros;
    r1(1:trialNumber) = zeros;
    acc(1:trialNumber,1:2) = zeros;
    rt(1:trialNumber,1:2) = zeros;
    Keyresponse(1:trialNumber,1:2) = zeros;
    
    %% instructions
    
    KbName('UnifyKeyNames');
%     crowdingInstructions(blockNum, w)
    Screen('FillRect', w, grey);
    Screen('TextSize', w, []);

    instText = [];

    %Practice Instructions
    if blockNum == 0 
        instText = ['Fixate on the center cross during the entire trial.\n'...
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
                    'You will hear a neutral beep when your response has been recorded.\n\n'...
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

    activeKeys = [KbName('LeftArrow') KbName('RightArrow') KbName('q')];
    RestrictKeysForKbCheck(activeKeys);
    
    %% main experiment
    
    while flag == 0 && trials < trialNumber % flag = 1 means that we've hit the upper limit
        trials = trials + 1;
        WhichStair = randi(hemifield*2); % which staircase to use
        stdir = stairdir(WhichStair);
        stori = ori(WhichStair);
        trial(WhichStair) = trial(WhichStair) + 1;  % count trials on this staircase
        
        priorityLevel = MaxPriority(w); % grab high priority to make generating movie as fast as possible
        Priority(priorityLevel);
        
        % show fixation cross
        Screen('FillRect', w, grey);
        % Draw fixation to indicate the start of the trial
        Screen('FillOval', w, white,FIXATION_POSITION,10);
        Screen('Flip',w);
        WaitSecs(markerWait(trials));
        
        % draw stimulus display
        % position the stimuli 
        dstRects(:,1) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials), yCen); % position the bars; target
        dstRects(:,2) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2)); 
            % above target
        dstRects(:,3) = CenterRectOnPoint(texrect, xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2)); 
            % below target
            
        % position the mask
        upperBound = dstRects(2,2);
        lowerBound = dstRects(4,3);
        leftBound = dstRects(1,3);
        rightBound = dstRects(3,2);
        nRow = 4;
        nCol = 4;  
        point = gridGen(nRow, nCol, upperBound, lowerBound, leftBound, rightBound);
        
        % pick orientations
        r1(trials) = oriIndex(trials)*ori(WhichStair);
        rotAngles = [r1(trials) -45 -45]; % optimal for threshold elevation
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles);
        % Screen('DrawText', w, sprintf('%g, stair = %d, correct = %d',r1(trials),WhichStair,stairCorrect(WhichStair)), ...
        %     30, 30, [255, 255, 255]);
        Screen('FillOval', w,white,FIXATION_POSITION,10);
        Screen('Flip',w);
        stimulus_onset_time(trials) = tic; % Mark the time when the display went up
        WaitSecs(presentationDur);
        
        % draw mask display
        for i = 1:length(point)
            maskDstRects(:,i) = CenterRectOnPoint(maskTexrect, point(i,1)+randi([-30 30]), point(i,2)+randi([-30 30]));
        end
        maskRotAngles = randi(360,1,length(point));
        Screen('FillOval', w,white,FIXATION_POSITION,10);
        Screen('DrawTextures', w, maskBarTex, [], maskDstRects, maskRotAngles);
        Screen('Flip',w);
        WaitSecs(maskDur);
        
        Screen('FillOval', w, white,FIXATION_POSITION,10);
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
                elseif responseKey(KbName('RightArrow'))
                    rspRatio(2) = rspRatio(2) + 1;
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
%         if blockNum
%             Beeper(750);
%         end
        if acc(trial(WhichStair),WhichStair) % IF CORRECT
            Beeper(1000); %play high beep for correct answer
            stairCorrect(WhichStair) = stairCorrect(WhichStair) + 1;
            if stairCorrect(WhichStair) == 3   % time to make stimulus harder?
                stairCorrect(WhichStair) = 0;  % zero counter
                
                % STAIRCASE DIRECTION FLAGS: 0 = down (getting harder), 1 = up (getting easier)
                if stairdir(WhichStair) == 1 % if stair direction is currently UP (= 1) this is a reversal
                    stairdir(WhichStair) = 0; % stair direction changes to DOWN
                    nReverse(WhichStair) = nReverse(WhichStair) + 1; % count reversal
                    stimulusReversal(WhichStair, nReverse(WhichStair)) = ori(WhichStair); % record stimulus value (stairStep) 
                                                                                            % at reversal
                end
                %change stimulus to make stimulus harder
                if ori(WhichStair) <= 5
                    ori(WhichStair) = ori(WhichStair) - delta;
                else
                    ori(WhichStair) = ori(WhichStair) - 2*delta;
                end
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
            if ori(WhichStair) <= 5
                ori(WhichStair)=ori(WhichStair) + delta;
            else
                ori(WhichStair) = ori(WhichStair) + 2*delta;
            end
            if ori(WhichStair) > maxori
                ori(WhichStair) = maxori;
            end % max of 45
        end
        
        ListenChar(0); %disables keyboard and flushes.
        ListenChar(2); %enables keyboard, no output to command window
        
        % test whether to end experiment
        if (nReverse(1) > 9 &&  nReverse(2) > 9) || trial(1) > 100 || trial(2) > 100
            flag = 1;
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
            stdir,stairCorrect(WhichStair),nReverse(WhichStair),hemiIndex(trials)],'-append', 'roffset', [],'delimiter', '\t');
        save(filenameTxt);
        save(filenameMatAll);
        save(filenameMat,'trials','trial','r1','acc','nReverse','stimulusReversal','rspRatio','hemiIndex');
        
        plot(stimulusReversal(1,1:nReverse(1)));hold on;
        plot(stimulusReversal(2,1:nReverse(2)));hold on;
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
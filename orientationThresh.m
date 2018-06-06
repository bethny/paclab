
%Makes a double staircase for Gabor orientation discrimination (one starts 
%at 0 degree difference, the other 10 degree) and after 10 reversals, 
%should show us at what point signal dots matter to the user. 
%3 down 1 up double staircase, estimate accuracy .792, d' = 1.634
% written by Jianfei, Fall 2015

atLab = 1;
directories = {'~/code','~/Bethany/paclab'};
addpath(genpath(directories{atLab+1}));

try
    VIEWING_DISTANCE_CM = 52;
    MONITOR_WIDTH_CM = 44;
    
%% subject information %%%%%%%%%%%%%%%%%%%    
    KbName('UnifyKeyNames')
    subjNum = input('\n Enter subject number: ');
    %Create a directory for this subject's data
    oripath=pwd;
    pathdata=strcat(pwd,filesep,'Subject_folders',filesep,'grasping_',subjNum,filesep);
    mkdir(pathdata);
    cd(pathdata);
    pathdata = pwd;
    thresholdFile = strcat(pathdata,filesep,subjNum,'_threshold.txt');
    cd(oripath);
    
%% open window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WhichScreen = max(Screen('Screens'));
    white = WhiteIndex(WhichScreen);
    grey = GrayIndex(WhichScreen);
    PsychDebugWindowConfiguration(0,0.75)
    
    [w, winRect] = Screen('OpenWindow',WhichScreen,128);
%     MyCLUT = load('gammaTable1.mat');
%     Screen('LoadNormalizedGammaTable', w, MyCLUT.gammaTable1*[1 1 1]);
    [xCen, yCen] = RectCenter(winRect);
    [swidth, sheight]=Screen('WindowSize', WhichScreen);
    screenCenter = [xCen yCen];
    screenInfo.center=screenCenter;
    degPerPixel = atan((MONITOR_WIDTH_CM/2)/VIEWING_DISTANCE_CM) * (180/pi) * (2 / swidth);
    ppd = 1 / degPerPixel;
    screenInfo.PPD = ppd;
    pixelPerCM = swidth / MONITOR_WIDTH_CM;
    CMperDeg=ppd/pixelPerCM;
    HideCursor;
    secPerFrame = Screen('GetFlipInterval',w);
                        
%% staircase value                      
    trials = 0; % count # trials overall?
    stairdir = [1 0];  % staircase 1 starts up (1) & 2 starts down (0) direction
    nReverse = [0 0];  % counts the number of reversals each staircase
    stairCorrect = [0 0];  % counts # correct in a row each staircase
    trial = [0 0]; % zero trial counters 2 staircases
    stimulusReversal(2,20) = zeros; % matrix with stimulus setting at staircase reversals 
    flag = 0;
    maxori = 10;
    minori = 0.5; % vertical
    ori = [minori maxori];    % 2 starting points
    delta = 1;          % how much to change signal at reversals\
    presentationDur = 1;
    
%% initial value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fixRad = 0.3;  %radius of fixation spot
    barWid = 0.4; % from bulakowski
    barLen = 4; % from bulakowski
    ecc = 10;
    tfDist = 4.3; % distance between target + flanker
    trialNumber = 200;
%     time1 = 0.13;  % duration of the 1st grating
%     time2 = 0.1;   % duration of ISI
    markerWaitList = [0.75, 1, 1.25];
    mn = 3;     % the number of markerWait;
    
    fixRadPx = round(fixRad*ppd);
    FIXATION_POSITION = [xCen - fixRadPx, yCen - fixRadPx, xCen + fixRadPx, yCen + fixRadPx];
    eccPx = round(ecc*ppd); 
    tfDistPx = round(tfDist*ppd);
    barLenPx = round(barLen*ppd);
    barWidPx = round(barWid*ppd);  

    markerWaitIndex = repmat(markerWaitList,1,floor((trialNumber+mn)./mn));
    n=randperm(trialNumber);
    markerWait= markerWaitIndex(1,n);
    
    barVert = ones(barLenPx, barWidPx); % assign pixel values for bar
    barTexVert = Screen('MakeTexture', w, barVert); % convert into texture
    texrect = Screen('Rect', barTexVert); % convert into rect
    dstRects(:,1) = CenterRectOnPoint(texrect, xCen + eccPx, yCen); % position the bars; target
    dstRects(:,2) = CenterRectOnPoint(texrect, xCen + eccPx + tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2)); % above target
    dstRects(:,3) = CenterRectOnPoint(texrect, xCen + eccPx - tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2)); % below target

    orientList = [-1 1]; % left vs right tilt
    orientIndex = repmat(orientList,1,trialNumber./2);
    n = randperm(trialNumber);
    oriIndex = orientIndex(1,n);
    
    %initialize stuff for feedback 
    FirstShow(1:trialNumber)=zeros;
    stimulus_onset_time(1:trialNumber)=zeros;
    r1(1:trialNumber)=zeros;
    acc(1:trialNumber,1:2)=zeros;
    rt(1:trialNumber,1:2)=zeros;
    Keyresponse(1:trialNumber,1:2)=zeros;
    
%% instruction

    % tell subjects what to look for
    % wait for keypress to continue
    KbName('UnifyKeyNames');
    Screen('DrawText', w, 'Press left arrow if the center bar is tilted to the left (counterclockwise)', 350, 400, [255, 255, 255]); 
    Screen('DrawText', w, 'Press right arrow if the center bar is tilted to the right (clockwise)',350, 500, [255, 255, 255]);
    Screen('Flip', w);
    KbWait;
    while KbCheck; end
    
    activeKeys = [KbName('LeftArrow') KbName('RightArrow') KbName('q')];
    RestrictKeysForKbCheck(activeKeys);

%% main experiment   

    while flag == 0 % flag = 1 means that we've hit the 45� upper limit
        trials = trials + 1;
        WhichStair = randi(2); % which staircase to use
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
        
%         fixShow = Screen('Flip',w);      
%         % Loop to ensure they stay at the marker for however long "markerWait" is
%         markerTime = tic;
%         while (toc(markerTime) < markerWait(1,trials))     
%         end     
        
        r1(trials) = oriIndex(trials)*ori(WhichStair);
        
        % draw stimulus display
        rotAngles = [r1(trials) -45 -45]; % optimal for threshold elevation
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles);
        Screen('DrawText', w, sprintf('%g, stair = %d, correct = %d',r1(trials),WhichStair,stairCorrect(WhichStair)), ...
            30, 30, [255, 255, 255]);
        Screen('FillOval', w,white,FIXATION_POSITION,10);
        
        Screen('Flip',w);
        WaitSecs(presentationDur);
        
%         FirstShow(trials) = Screen('Flip',w,fixShow+markerWait(trials));
        stimulus_onset_time(trials) = tic; % Mark the time when the display went up
        
        Screen('FillOval', w, white,FIXATION_POSITION,10);
        Screen('Flip',w);
%         ISITime(trials) = Screen('Flip',w,FirstShow(trials)+presentationDur);
       
        keypressed = 0;
        while (keypressed==0) 
            [keyIsDown, secs, responseKey] = KbCheck; 

            if keyIsDown
                if responseKey(KbName('ESCAPE'))
                    Screen('CloseAll');
                    ShowCursor;
                end

                if oriIndex(trials) == 1 && responseKey(KbName('RightArrow')) || oriIndex(trials) == -1 && ...
                        responseKey(KbName('LeftArrow'))
                    acc(trial(WhichStair),WhichStair) = 1;
                    rt(trial(WhichStair),WhichStair) = (secs-stimulus_onset_time(trials));
                    Keyresponse(trial(WhichStair),WhichStair)=find(responseKey); 
                else
                    acc(trial(WhichStair),WhichStair) = 0;
                    rt(trial(WhichStair),WhichStair) = (secs-stimulus_onset_time(trials));
                    Keyresponse(trial(WhichStair),WhichStair)=find(responseKey); 
                end
                KbReleaseWait;              
                keypressed=1;
            end
        end    
        
       %staircase stuff        
       if acc(trial(WhichStair),WhichStair) % IF CORRECT
           Beeper(1000); %play high beep for correct answer
           stairCorrect(WhichStair) = stairCorrect(WhichStair) + 1;
           if stairCorrect(WhichStair) == 3   % time to make stimulus harder?
               stairCorrect(WhichStair) = 0;  % zero counter
               
               % STAIRCASE DIRECTION FLAGS: 0 = down (getting harder), 1 = up (getting easier)
               if stairdir(WhichStair) == 1 % if stair direction is currently UP (= 1) this is a reversal
                   stairdir(WhichStair) = 0; % stair direction changes to DOWN
                   nReverse(WhichStair) = nReverse(WhichStair) + 1; % count reversal
                   stimulusReversal(WhichStair, nReverse(WhichStair)) = ori(WhichStair); % record stimulus value (stairStep) at reversal
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
       dlmwrite(thresholdFile,[WhichStair,trial(WhichStair),stori,r1(trials),acc(trial(WhichStair),WhichStair),...
            stdir,stairCorrect(WhichStair),nReverse(WhichStair)],'-append', 'roffset', [],'delimiter', '\t');
       
       ListenChar(0); %disables keyboard and flushes.
       ListenChar(2); %enables keyboard, no output to command window
       
       % test whether to end experiment
       if (nReverse(1) > 9 &&  nReverse(2) > 9) || trial(1) > 100 || trial(2) > 100 
           flag = 1; 
       end
       save('threshold.mat')
    end
    
    Screen ('CloseAll');
	ShowCursor;
    %Priority(0);

%--------- DONE --------------------
% rudimentary data analysis
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

    pathdata=strcat(pwd,filesep,'Subject_folders',filesep,'grasping_',subjNum,filesep);
    mkdir(pathdata);
    cd(pathdata);
    save('threshold.mat')
    plot(stimulusReversal(1,1:nReverse(1)));hold on;
    plot(stimulusReversal(2,1:nReverse(2)));hold on;
    disp 
catch psychlasterror
    Screen ('CloseAll');
    close all;
    ShowCursor;
    disp(psychlasterror.message);
    disp(psychlasterror.stack);
    %psychrethrow (psychlasterror);

end   % try catch
% end   % function


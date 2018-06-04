                                                                                                                                                                                 function  discrimination_threshold_single
%Makes a double staircase for Gabor orientation discrimination (one starts 
%at 0 degree difference, the other 10 degree) and after 10 reversals, 
%should show us at what point signal dots matter to the user. 
%3 down 1 up double staircase, estimate accuracy .792, d' = 1.634
% written by Jianfei, Fall 2015

atLab = 1;
directories = {'~/code','~/Bethany/paclab'};
addpath(genpath(directories{atLab+1}));

try
    %rand('twister', sum(100*clock)); % seed RNG
    VIEWING_DISTANCE_CM = 52;%57;
    MONITOR_WIDTH_CM = 44;%35;
    
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
    PsychDebugWindowConfiguration(0,0.5)
    
    PsychDefaultSetup(2);
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
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
    nReverse = [0 0 0];  % counts the number of reversals each staircase
    stairCorrect = [0 0 0];  % counts # correct in a row each staircase
    trial = [0 0 0]; % zero trial counters 2 staircases
    stimulusReversal(2,20) = zeros; % matrix with stimulus setting at staircase reversals 
    flag = 0;
    maxori = 45;
    minori = 0; % vertical
    ori = [minori maxori];    % orientation difference for 2 staircases / is this 2 starting points? unclear
    delta = 1;          % how much to change signal & noise dot numbers at reversals?
    
%% initial value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fixRad = 0.7/2;  %radius of fixation spot
    barWid = 0.4; % from bulakowski
    barLen = 4; % from bulakowski
    ecc = 9;
    tfDist = 4.3; % distance between target + flanker
    startOri = 45/2;
    orientation2 = [45 135];
    orin2 = 2;     % the number of orientation2
    trialNumber = 350;
    time1 = 0.13;  % duration of the 1st grating
    time2 = 0.1;   % duration of ISI
    markerWaitList = [0.75, 1, 1.25];
    mn=3;     % the number of markerWait;
    pn=2;     % the number of position

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

    orientList2=1:orin2;
    orientIndex2 = repmat(orientList2,1,trialNumber./orin2);
    n=randperm(trialNumber);
    oriIndex2 = orientIndex2(1,n);

    orientList1=[-1 1];
    orientIndex1 = repmat(orientList1,1,trialNumber./2);
    n=randperm(trialNumber);
    oriIndex1 = orientIndex1(1,n);

    conditionList1=[1 1 2 2 2];
    conditionIndex1 = repmat(conditionList1,1,trialNumber./5);
    n=randperm(trialNumber);
    conIndex = conditionIndex1(1,n);
    
    %initialize stuff for feedback 
    %high tone for correct and low tone for incorrect
    %         based on Erika's code (& Max's code) in partial Italian
    fs = 44100;             %no idea what this is for
    durata = .4;            %duration of tone?
    frequenza_high = 1000;  %high tone 1000 Hz?
    frequenza_low = 500;    %low tone 500 Hz?
    intensita = .5;         %amplitude?
    segnale_high = sin((1/fs:1/fs:durata) * 2 * pi * frequenza_high) * intensita;
    segnale_low = sin((1/fs:1/fs:durata) * 2 * pi * frequenza_low) * intensita; 
    
    FirstShow(1:trialNumber)=zeros;
    SecondShow(1:trialNumber)=zeros;
    stimulus_onset_time(1:trialNumber)=zeros;
    r1(1:trialNumber)=zeros;
    r2(1:trialNumber)=zeros;
    r3(1:trialNumber)=zeros;
    acc(1:trialNumber,1:3)=zeros;
    rt(1:trialNumber,1:3)=zeros;
    Keyresponse(1:trialNumber,1:3)=zeros;
%% instruction

    % tell subjects what to look for
    % wait for keypress to continue
    KbName('UnifyKeyNames');
    activeKeys = [KbName('LeftArrow') KbName('RightArrow') KbName('q')];
    RestrictKeysForKbCheck(activeKeys);
    srcRect=[];
    Screen('DrawText', w, 'Press left arrow if the center bar is tilted to the left (counterclockwise)', 350, 400, [255, 255, 255]); 
    Screen('DrawText', w, 'Press right arrow if the center bar is tilted to the right (clockwise)',350, 500, [255, 255, 255]);
    Screen('Flip', w);
    KbWait;
    while KbCheck; end;

%% main experiment   

    while flag == 0 % flag = 1 means that we've hit the 45º upper limit
        trials = trials + 1;
        if conIndex(trials)==2 % if conIdx = 2, use up or down staircases (WhichStair = 1 || 2)
            WhichStair = randi(2); % which staircase to use
            stdir=stairdir(WhichStair);
            stori=ori(WhichStair);
        else
            WhichStair = 3; % which staircase to use
            stori=0;
            stdir=-1;
        end
        trial(WhichStair) = trial(WhichStair) + 1;  % count trials on this staircase
                
        priorityLevel=MaxPriority(w); % grab high priority to make generating movie as fast as possible
        Priority(priorityLevel);

        % show fixation cross
        Screen('FillRect', w, grey);
        % Draw fixation to indicate the start of the trial       
        Screen('FillOval', w, white,FIXATION_POSITION,10);     
        Screen('Flip',w);
        % Loop to ensure they stay at the marker for however long "markerWait" is
        markerTime = tic;
        while (toc(markerTime) < markerWait(1,trials))
     
        end     
        if conIndex(trials)==2            
            r1(trials)=oriIndex1(trials)*ori(WhichStair); % oriIndex1 = left or right tilt; ori = 
            r2(trials)=orientation2(oriIndex2(trials));
            r3(trials)=r1(trials)+r2(trials);
        else
            r1(trials)=0;
            r2(trials)=orientation2(oriIndex2(trials));
            r3(trials)=r1(trials)+r2(trials);
        end
        rl=barWidPx./sqrt(2);
        rv=VgaborDimPix./sqrt(8);
        if r2(trials)==45
            lp1=[pos(dstIndex(trials),1)-rl,pos(dstIndex(trials),2)-rl,pos(dstIndex(trials),1)-rv,pos(dstIndex(trials),2)-rv];
            lp2=[pos(dstIndex(trials),1)+rv,pos(dstIndex(trials),2)+rv,pos(dstIndex(trials),1)+rl,pos(dstIndex(trials),2)+rl];
        else
            lp1=[pos(dstIndex(trials),1)+rv,pos(dstIndex(trials),2)-rv,pos(dstIndex(trials),1)+rl,pos(dstIndex(trials),2)-rl];
            lp2=[pos(dstIndex(trials),1)-rl,pos(dstIndex(trials),2)+rl,pos(dstIndex(trials),1)-rv,pos(dstIndex(trials),2)+rv];
        end
            % draw stimulus display
        rotAngles = [ori -45 -45]; % optimal for threshold elevation
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles);
%         Screen('LineStipple',w,1,2);
%         Screen('DrawLine', w,[],lp1(1),lp1(2),lp1(3),lp1(4));
%         Screen('DrawLine', w,[],lp2(1),lp2(2),lp2(3),lp2(4));
        
        Screen('FillOval', w,white,FIXATION_POSITION,10);
        

        FirstShow(trials)=Screen('Flip',w);

        Screen('FillOval', w, white,FIXATION_POSITION,10);     
        ISITime(trials)=Screen('Flip',w,FirstShow(trials)+time1);
      
        % Mark the time when the display went up
        stimulus_onset_time(trials) = tic;
        
        keypressed=0;
        while (keypressed==0) 
            [keyIsDown, secs, responseKey]=KbCheck; 

            if keyIsDown
                if responseKey(KbName('ESCAPE'))
                    Screen('CloseAll');
                    ShowCursor;
                end
                if (conIndex(trials)==2) && (find(responseKey)==key2)
                    acc(trial(WhichStair),WhichStair)=1;
                    rt(trial(WhichStair),WhichStair) = (secs-stimulus_onset_time(trials));
                    Keyresponse(trial(WhichStair),WhichStair)=find(responseKey); 
                elseif (conIndex(trials)==2) && (find(responseKey)==key1)
                    acc(trial(WhichStair),WhichStair)=0;
                    rt(trial(WhichStair),WhichStair) = (secs-stimulus_onset_time(trials));
                    Keyresponse(trial(WhichStair),WhichStair)=find(responseKey); 
                elseif (conIndex(trials)==1) && (find(responseKey)==key1)
                    acc(trial(3),3)=1;
                    rt(trial(3),3) = (secs-stimulus_onset_time(trials));
                    Keyresponse(trial(3),3)=find(responseKey); 
                elseif (conIndex(trials)==1) && (find(responseKey)==key2)
                    acc(trial(3),3)=0;
                    rt(trial(3),3) = (secs-stimulus_onset_time(trials));
                    Keyresponse(trial(3),3)=find(responseKey);
                end                        
                KbReleaseWait;              
                keypressed=1;
            end;

        end        

       %staircase stuff        
       if (conIndex(trials)==2)&& acc(trial(WhichStair),WhichStair)
            sound ( segnale_high, fs); %play high beep for correct answer
            stairCorrect(WhichStair) = stairCorrect(WhichStair) + 1;
            if stairCorrect(WhichStair) == 3   % time to make stimulus harder?
                stairCorrect(WhichStair) = 0;  % zero counter

    % STAIRCASE DIRECTION FLAGS: 0 = down (getting harder), 1 = up (getting easier) 

                % if stair direction is currently UP (= 1) this is a reversal
                if stairdir(WhichStair) == 1

                    stairdir(WhichStair) = 0; % stair direction changes to DOWN
                    nReverse(WhichStair) = nReverse(WhichStair) + 1; % count reversal 
                     % record stimulus value (stairStep) at reversal
                    stimulusReversal(WhichStair, nReverse(WhichStair)) = ori(WhichStair); 
                end 
                %change stimulus/noise dots (make there be more noise dots,
                %fewer signal dots (+/- delta) to make stimulus harder
                if (ori(WhichStair)<=10)
                    ori(WhichStair)=ori(WhichStair)-delta;
                else                       
                ori(WhichStair) = ori(WhichStair) - 2*delta;
                end
                if ori(WhichStair) < minori, ori(WhichStair) = minori; end %there can't be negative dots
            end 


       elseif (conIndex(trials)==2) && acc(trial(WhichStair),WhichStair)==0 %incorrect response
           sound ( segnale_low, fs); %play low beep for incorrect answer
           stairCorrect(WhichStair) = 0;
           if stairdir(WhichStair) == 0                             
              stairdir(WhichStair) = 1;   % change direction to UP
              nReverse(WhichStair) = nReverse(WhichStair) + 1; % count reversal
              stimulusReversal(WhichStair, nReverse(WhichStair)) = ori(WhichStair); % record stimulus @ this reversal
           end
             %change stimulus/noise dots (make there be fewer noise dots,
            %more signal dots (+/- 10) to make stimulus easier
           if (ori(WhichStair)<=10)
                ori(WhichStair)=ori(WhichStair)+delta;
           else                       
            ori(WhichStair) = ori(WhichStair) + 2*delta;
           end           
           if ori(WhichStair) > maxori, ori(WhichStair) = maxori; end % max of 45 
       elseif (conIndex(trials)==1)&& acc(trial(3),3)
           sound ( segnale_high, fs);
           stairCorrect(3) = -1;
           nReverse(3)=-1;
       elseif (conIndex(trials)==1)&& acc(trial(3),3)==0
           sound ( segnale_low, fs)
           stairCorrect(3) = -1;
           nReverse(3)=-1;
       end;
        dlmwrite(thresholdFile,[WhichStair,trial(WhichStair),stori,r2(trials),acc(trial(WhichStair),WhichStair),stdir,stairCorrect(WhichStair),nReverse(WhichStair),dstIndex(trials)],'-append', 'roffset', [],'delimiter', '\t');
        
         ListenChar(0); %disables keyboard and flushes. 
         ListenChar(2); %enables keyboard, no output to command window
        
         % test whether to end experiment
        if (nReverse(1) > 9 &&  nReverse(2) > 9) || trial(1) > 100 || trial(2) > 100; flag = 1; end
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

    % whos   
     %close all;
     %clear all;
     %threshold=[stair1mean,stair2mean,StandardDev1,StandardDev2;(stair1mean+stair2mean)./2,(StandardDev1+StandardDev2)./2];
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
end   % function


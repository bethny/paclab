function [output] = gaborDisplay

%%  TO-DO LIST
%   1. Add masking after ~100ms

Screen('Preference', 'SkipSyncTests', 1); % only for testing purposes

atLab = 1;
directories = {'~/code','~/Bethany/paclab'};
addpath(genpath(directories{atLab+1}));

fixationDur = .75; % duration of fixation cross only
presentationDur = .2; % duration gabors are visible

%% SCREEN PARAMETERS
viewDist = 52;
monitorWidth = 44;
screenid = max(Screen('Screens'));
[swidth, ~] = Screen('WindowSize', screenid);
degPerPixel = atan((monitorWidth/2)/viewDist) * (180/pi) * (2/swidth);
ppd = 1/degPerPixel;

%% USER INPUT
subj = input('Subject name: ','s');
block = input('Block number (0 = practice, 1-x = exp): ');
dataDir = {'~/code/pac/Data/Psychophysics','~/Bethany/paclab/Data/Psychophysics'};
fileDir = sprintf('%s/%s',dataDir{atLab+1},subj);
if ~exist(fileDir)
    mkdir(fileDir)
end

%% GABOR PATCH PARAMETERS

% LIVNE & SAGI PARAMETERS
ecc = 10;
radialDist = 4;

% JIANFEI'S PARAMETERS
contrast = 1;
gsize = 5;      %size of the grating
vsize = 2.5;      % LIVNE & SAGI / size of the visible grating (edge size)
% reachsize = 3;  %size of the reachable area (the threshold for subject grasping)
% degree = 8.5;    %distance between center and grating
sf = 1;       %spatial frequency

eccPix = ceil(ecc*ppd);
radialDistPix = ceil(radialDist*ppd);
gaborDimPix = round(gsize*ppd);
VgaborDimPix = round(vsize*ppd);
sigma = gaborDimPix/8;
numCycles = gsize*sf;
freq = numCycles / gaborDimPix;     %spatial frequency of the grating
freqmask = 0.00001;

%% LEVELS
% define all conditions
crowding = 0:2;
tilt = [-40 -30 -20 -10 10 20 30 40];
hemifield = [-eccPix eccPix];

nRepeat = 1; % # of times to repeat each cnd

[p,q,r] = meshgrid(crowding, tilt, hemifield);
cnd = [p(:) q(:) r(:)];

% randomize & replicate conditions for full list of stimuli
allTrials = repmat(cnd,[nRepeat,1]);
shuffledTrials = allTrials(randperm(size(allTrials,1)),:); 

if strcmp(subj,'test')
    nTrials = 6;
elseif block == 0
    nTrials = 20;
else
    nTrials = length(shuffledTrials);
end

%% OPENING PTB
PsychDebugWindowConfiguration(0,1)
PsychDefaultSetup(2);

PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
[win, winRect] = PsychImaging('OpenWindow', screenid, 0.5);
    
KbQueueCreate; %creates cue using defaults
KbQueueStart;  %starts the cue

%% TRIALS
KbName('UnifyKeyNames');
activeKeys = [KbName('LeftArrow') KbName('RightArrow') KbName('q')];
t2wait = 2; % maximum time to wait for response
RestrictKeysForKbCheck(activeKeys);
ListenChar(2);
white = WhiteIndex(win);
black = BlackIndex(win);
grey = GrayIndex(win, 0.5);

HideCursor
Screen(win, 'FillRect', grey);
Screen('Flip', win); 
Screen('TextSize', win, 20);

crowdingInstructions(block, win)

for i = 1:nTrials
    clear texrect
    clear texrectMask
    clear inrect
    clear inrectMask
    clear dstRects
    clear dstRectsMask
    
    crowding = shuffledTrials(i,1);
    tilt = shuffledTrials(i,2);
    eccPix = shuffledTrials(i,3);

    if ~crowding
        nGabors = 1;
    else
        nGabors = 7;
    end

    [w, h] = RectSize(winRect);
    ifi = Screen('GetFlipInterval', win); % how necessary is this? no animations
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);

    m = grey*contrast*makegabor(gaborDimPix,freq,sigma,0,VgaborDimPix); %make gabor
    gabortex = Screen('MakeTexture', win, m,[],[],2);
    texrect = Screen('Rect', gabortex);
    texrectMask = [0 0 VgaborDimPix VgaborDimPix];
    inrect = repmat(texrect', 1, nGabors);
    inrectMask = repmat(texrectMask',1,nGabors);
     
    % STIMULUS PARAMETERS
    dstRects = zeros(4, nGabors); 
    scale = ones(1,nGabors);
    refPtX = w/2 + eccPix;
    refPtY = h/2;

    % TARGET GABOR LOCATION
    dstRects(:,1) = CenterRectOnPoint(texrect, refPtX, refPtY)';
    dstRectsMask(:,1) = CenterRectOnPoint(texrectMask, refPtX, refPtY)';

    % FLANKER GABORS LOCATION
    if crowding        
        dstRects(:,2) = CenterRectOnPoint(texrect, refPtX + radialDistPix, refPtY)'; % rightmost
        dstRects(:,3) = CenterRectOnPoint(texrect, refPtX + radialDistPix*cos(deg2rad(60)), refPtY - radialDistPix*sin(deg2rad(60)))';
        dstRects(:,4) = CenterRectOnPoint(texrect, refPtX - radialDistPix*cos(deg2rad(60)), refPtY - radialDistPix*sin(deg2rad(60)))';
        dstRects(:,5) = CenterRectOnPoint(texrect, refPtX - radialDistPix, refPtY)'; % leftmost
        dstRects(:,6) = CenterRectOnPoint(texrect, refPtX - radialDistPix*cos(deg2rad(60)), refPtY + radialDistPix*sin(deg2rad(60)))';
        dstRects(:,7) = CenterRectOnPoint(texrect, refPtX + radialDistPix*cos(deg2rad(60)), refPtY + radialDistPix*sin(deg2rad(60)))';
        dstRectsMask(:,2) = CenterRectOnPoint(texrectMask, refPtX + radialDistPix, refPtY)'; % rightmost
        dstRectsMask(:,3) = CenterRectOnPoint(texrectMask, refPtX + radialDistPix*cos(deg2rad(60)), refPtY - radialDistPix*sin(deg2rad(60)))';
        dstRectsMask(:,4) = CenterRectOnPoint(texrectMask, refPtX - radialDistPix*cos(deg2rad(60)), refPtY - radialDistPix*sin(deg2rad(60)))';
        dstRectsMask(:,5) = CenterRectOnPoint(texrectMask, refPtX - radialDistPix, refPtY)'; % leftmost
        dstRectsMask(:,6) = CenterRectOnPoint(texrectMask, refPtX - radialDistPix*cos(deg2rad(60)), refPtY + radialDistPix*sin(deg2rad(60)))';
        dstRectsMask(:,7) = CenterRectOnPoint(texrectMask, refPtX + radialDistPix*cos(deg2rad(60)), refPtY + radialDistPix*sin(deg2rad(60)))';
    end

    % ROTATION ANGLES
    % 0 is vertical ; angles go CW
    flankerAngle = [0 120 240]; 
    if ~crowding
        rotAngles = tilt;
    elseif crowding == 1
        rotAngles = [tilt repmat(flankerAngle,[1,2])];
    else
        rotAngles = [tilt repmat(flankerAngle,[1,2])+90];
    end
    
    % BLANK SCREEN ITI 
    trialITI = rand+.5;
    WaitSecs(trialITI) %jitters prestim interval between .5 and 1.5 seconds; inter-trial interval  

    % DRAW FIXATION CROSS OR TRIAL NUMBER
    fixCrossDim = 0.5;
    fixCrossDimPix = fixCrossDim*ppd; % size
    xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
    yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
    allCoords = [xCoords; yCoords];
    lineWidthPix = 1; % linewidth
    
    if ~mod(nTrials-i,10)
        DrawFormattedText(win, sprintf('%d',nTrials-i), 'Center', 'Center', [255 255 255]);
    else
        
        Screen('DrawLines', win, allCoords,lineWidthPix, WhiteIndex(screenid), [w/2 h/2], 2);
    end
    
    Screen('Flip',win); % PRESENT FIXATION ONLY
    WaitSecs(fixationDur)
      
    % DRAW GABOR PATCHES + FIXATION 
    Screen('DrawLines', win, allCoords,lineWidthPix, WhiteIndex(screenid), [w/2 h/2], 2);
    Screen('DrawTextures', win, gabortex, [], dstRects, rotAngles, [], 0.5, [], [], 0);
    Screen('DrawingFinished', win);
    
    % CENTER GABORS
    [x, y] = RectCenterd(dstRects);
    if ~crowding
        dstRects = CenterRectOnPointd((inrect .* repmat(scale,4,1))', x, y); 
    else
        dstRects = CenterRectOnPointd(inrect .* repmat(scale,4,1), x, y);
    end

    vbl = Screen('Flip', win); % PRESENT GABORS
    WaitSecs(presentationDur)
    
    % FIXATION & MASK
    Screen('DrawLines', win, allCoords,lineWidthPix, WhiteIndex(screenid), [w/2 h/2], 2);
    Screen('FillOval', win, [], dstRectsMask);
    Screen('Flip', win); % PRESENT MASK
    
    trialStart = tic;   
    timedout = 1;
    
    while toc(trialStart) < t2wait && timedout
        [keyIsDown, keyTime, keyCode] = KbCheck;
        if keyIsDown && ~strcmp(KbName(keyCode),'q')
            Beeper;
            timedout = 0;
            timeElapsed = toc(trialStart);
        elseif keyIsDown && strcmp(KbName(keyCode),'q')
            ShowCursor;
            ListenChar(0);
            Screen('CloseAll');
            return;
        end
    end
    
    if ~timedout
        rsp.RT(i) = timeElapsed;
        rsp.keyName{i} = KbName(keyCode);
    else
        Beeper(900);
        rsp.RT(i) = t2wait;
        rsp.keyName{i} = 'none';
    end

    vbl = Screen('Flip', win, vbl + 0.5 * ifi); 
    
end
ListenChar(1)
ShowCursor();
sca

output.block = block;
output.subj = subj;
output.rsp = rsp;
output.cnd = shuffledTrials;

if block ~= 0 % RE-ENABLE LATER
    save(sprintf('%s/%s_%d.mat',fileDir,subj,block),'output');
    if ~exist(sprintf('%s/Raw',fileDir))
        mkdir(sprintf('%s/Raw',fileDir))
    end
    save(sprintf('%s/Raw/%s_%d_raw.mat',fileDir,subj,block));
end

return

%     
%     if ~rspGiven
%         Beeper(900);
%         timedout = 1; 
%     end
%     
%     
%     stimulus_onset_time = tic;    
%     response(i) = -1;
%     % Loop to look for reach movements
%     while (toc(stimulus_onset_time) < RTDeadline && ~exit)
% 
%         [keyIsDown, secs, keyCode]= KbCheck;
%           key = KbName(keyCode);
%             if ((keyIsDown)&& ((strcmp(key(1),'z')) || (strcmp(key(1),'x')) || (strcmp(key(1),'c')) || (strcmp(key(1),'v'))))
%                timeElapsed = toc(stimulus_onset_time);  
%                if strcmp(key(1),'z')
%                    response(i) = 1;
%                elseif  strcmp(key(1),'x')
%                    response(i) = 2;
%                elseif strcmp(key(1),'c')
%                    response(i) = 3;
%                elseif strcmp(key(1),'v')
%                    response(i) = 4;
%                end
%                 exit = 1;
%             elseif ((keyIsDown)&& strcmp(key(1),'q'))
%                         showcursor;
%                         ListenChar(0);
%                         Screen('CloseAll');
%                         return;
%             end
% 
%      
%     end
%     
%     %%%%%%% START %%%%%%%
%     
%     
%     tStart = GetSecs;
%     timedout = false;
% %     while ~timedout
% %         [keyIsDown, keyTime, keyCode] = KbCheck; 
% %         if keyIsDown
% %             Beeper;
% %             break; 
% %         end
% %         if (keyTime - tStart) > t2wait
% %             Beeper(900);
% %             timedout = true; 
% %         end
% %     end
%     
%     if ~timedout
%         rsp.RT(i) = keyTime - tStart;
%         rsp.keyName{i} = KbName(keyCode);
%     else
%         rsp.RT(i) = t2wait;
%         rsp.keyName{i} = 'none';
%     end
    %%%%%%% END %%%%%%%
function [output] = gaborDisplay

%%  TO-DO LIST
%   1. Compute visual angle based on viewing distance & size of screen
%   3. Look up best timing parameters. 
%   4. Reject keypresses during fixation-cross-only display.

addpath(genpath('~/code/pac'));

fixationDur = 1;

%% USER INPUT
subj = input('Subject name: ','s');
block = input('Block number (0 = practice, 1-x = exp): ');
fileDir = sprintf('~/code/pac/Data/Psychophysics/%s',subj);
if ~exist(fileDir)
    mkdir(fileDir)
end

if subj == 'test'
    nTrials = 6;
elseif block == 0
    nTrials = 20;
else
    nTrials = 144;
end

%% LEVELSz
crowding = [0:2];
tilt = [-40 -30 -20 -10 10 20 30 40];

[p,q] = meshgrid(crowding, tilt);
pairs = [p(:) q(:)];

allPairs = repmat(pairs,[6,1]);
shuffledPairs = allPairs(randperm(size(allPairs,1)),:);

%% GABOR PATCH PARAMETERS

% DIMENSIONS IN VISUAL DEGREES
gaborDiam = 2;

% Pixel per º/vis ang
ppd = 40.5;

% si = 50; 
si = ceil(ppd*(gaborDiam/2))+1; % half the desired gabor size in PIXELS

% Size of support in pixels, derived from si:
tw = 2*si;
th = 2*si;

phase = 0.5; % Phase of underlying sine grating in degrees:
sc = 18.0; % Spatial constant of the exponential "hull"
freq = 0.05;
contrast = 25;
aspectratio = 1.0; % width / height

PsychDefaultSetup(2);
screenid = max(Screen('Screens'));
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
grey = GrayIndex(win,0.5);
HideCursor
Screen(win, 'FillRect', grey);
Screen('Flip', win); 
Screen('TextSize', win, 20);

crowdingInstructions(block, win)

for i = 1:nTrials
    crowding = shuffledPairs(i,1);
    tilt = shuffledPairs(i,2);

    if ~crowding
        nGabors = 1;
    else
        nGabors = 7;
    end

    [w, h] = RectSize(winRect);
    ifi = Screen('GetFlipInterval', win); % how necessary is this? no animations
    Screen('BlendFunction', win, GL_ONE, GL_ONE);
    mypars = repmat([phase+180, freq, sc, contrast, aspectratio, 0, 0, 0]', 1, nGabors);
    gabortex = CreateProceduralGabor(win, tw, th, 1);
    texrect = Screen('Rect', gabortex);
    inrect = repmat(texrect', 1, nGabors);
    
    ecc = w/4; % eccentricity of target gabor

    % STIMULUS PARAMETERS
    dstRects = zeros(4, nGabors); 
    scale = ones(1,nGabors);
    refPtX = w/2 + ecc;
    refPtY = h/2;

    % TARGET GABOR LOCATION
    dstRects(:,1) = CenterRectOnPoint(texrect, refPtX, refPtY)';

    % FLANKER GABORS LOCATION
    if crowding
        radialDist = si*4;
        dstRects(:,2) = CenterRectOnPoint(texrect, refPtX + radialDist, refPtY)';
        dstRects(:,3) = CenterRectOnPoint(texrect, refPtX + radialDist*cos(deg2rad(60)), refPtY - radialDist*sin(deg2rad(60)))';
        dstRects(:,4) = CenterRectOnPoint(texrect, refPtX - radialDist*cos(deg2rad(60)), refPtY - radialDist*sin(deg2rad(60)))';
        dstRects(:,5) = CenterRectOnPoint(texrect, refPtX - radialDist, refPtY)';
        dstRects(:,6) = CenterRectOnPoint(texrect, refPtX - radialDist*cos(deg2rad(60)), refPtY + radialDist*sin(deg2rad(60)))';
        dstRects(:,7) = CenterRectOnPoint(texrect, refPtX + radialDist*cos(deg2rad(60)), refPtY + radialDist*sin(deg2rad(60)))';
    end

    % ROTATION ANGLES
    % 0 is vertical; angles go CCW
    flankerAngle = [0 60 120]; 
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

    % DRAW FIXATION CROSS
    fixCrossDimPix = 20; % size
    xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
    yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
    allCoords = [xCoords; yCoords];
    lineWidthPix = 1; % linewidth
    Screen('DrawLines', win, allCoords,lineWidthPix, WhiteIndex(screenid), [w/2 h/2], 2);
    Screen('Flip',win);
    WaitSecs(fixationDur)
      
    % DRAW GABOR PATCHES + FIXATION 
    Screen('DrawLines', win, allCoords,lineWidthPix, WhiteIndex(screenid), [w/2 h/2], 2);
    Screen('DrawTextures', win, gabortex, [], dstRects, rotAngles, [], [], [], [], kPsychDontDoRotation, mypars);
    Screen('DrawingFinished', win);
    
    % CENTER GABORS
    [x, y] = RectCenterd(dstRects);
    if ~crowding
        dstRects = CenterRectOnPointd((inrect .* repmat(scale,4,1))', x, y); 
    else
        dstRects = CenterRectOnPointd(inrect .* repmat(scale,4,1), x, y);
    end

    vbl = Screen('Flip', win);    
    trialStart = tic;   
    timedout = 1;
    
    while toc(trialStart) < t2wait && timedout
        [keyIsDown, keyTime, keyCode] = KbCheck;
        if keyIsDown
            Beeper;
            timedout = 0;
            timeElapsed = toc(trialStart);
        elseif keyIsDown && strcmp(key(1),'q')
            showcursor;
            ListenChar(0);
            Screen('CloseAll');
            return;
        end
    end
    
    if ~timedout
        rsp.RT(i) = timeElapsed;
        rsp.keyName{i} = KbName(keyCode);
    else
        rsp.RT(i) = t2wait;
        rsp.keyName{i} = 'none';
    end
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
    vbl = Screen('Flip', win, vbl + 0.5 * ifi); 
    
end
ListenChar(1)
ShowCursor();
sca

output.block = block;
output.subj = subj;
output.rsp = rsp;
output.cnd = shuffledPairs;

if block ~= 0
    save(sprintf('%s/%s_%d.mat',fileDir,subj,block),'output');
    if ~exist(sprintf('%s/Raw',fileDir))
        mkdir(sprintf('%s/Raw',fileDir))
    end
    save(sprintf('%s/Raw/%s_%d_raw.mat',fileDir,subj,block));
end

return
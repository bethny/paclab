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
    clear all
    close all
    
    VIEWING_DISTANCE_CM = 52;
    MONITOR_WIDTH_CM = 44;
    
    %% subject information %%%%%%%%%%%%%%%%%%%
    KbName('UnifyKeyNames')
    
    subjNum = input('\n Enter subject number: ');
    blockNum = input('\n Enter block number: ');
    resume = input('\n Resume a crashed experiment? (0 = no, 1 = yes): ');
    grasping = input('\n No-action (0) or grasp (1): ');
    task = {'noAction','grasp'};
   
    oripath = 'C:\Documents and Settings\js21\My Documents\MATLAB\Bethany';
    addpath(genpath(oripath));
    cd(oripath);
    pathdata = strcat(pwd,filesep,'Subject_folders',filesep,'S0',num2str(subjNum),filesep);
    
    if grasping
        calibrate = input('\n Calibrate? (0 = no, 1 = yes): ');
        calibFile = strcat(pathdata,num2str(subjNum),'_calibration.mat');  %% calibration filename
        if exist(calibFile, 'file') ~= 2 || calibrate
            [affine_alignment1, screen_y1, marker_pos1] = calibrate_sensor(1);
            [affine_alignment2, screen_y2, marker_pos2] = calibrate_sensor(2);
            save(calibFile,'affine_alignment1','screen_y1','marker_pos1','affine_alignment2','screen_y2','marker_pos2');
        else
            load(calibFile);
        end
    end
    
    if blockNum
        cd(pathdata);
        filenameMat = strcat(pathdata,sprintf('%dblock%d_%s',subjNum,blockNum,task{grasping+1}),'.mat');
        filenameMatAll = strcat(pathdata,sprintf('%dblock%d_%s',subjNum,blockNum,task{grasping+1}),'_all.mat');
    end
    
    if resume 
        load(filenameMatAll);
        fprintf(sprintf('Experiment resumed at trial %d\n',trials))
        resume = 1;
    end
    
    % get threshold info
    data = load(sprintf('%s%dblock1_threshold_all.mat',pathdata,subjNum));
    thresholds = data.thresholds;
    
    %% open window %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WhichScreen = max(Screen('Screens'));
    white = WhiteIndex(WhichScreen);
    grey = GrayIndex(WhichScreen);
    
    [w, winRect] = Screen('OpenWindow',WhichScreen,128);
    MyCLUT = load('C:\Documents and Settings\js21\My Documents\MATLAB\Bethany\gammaTable1.mat');
    Screen('LoadNormalizedGammaTable', w, MyCLUT.gammaTable1*[1 1 1]);

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

    %% Tracker Thresholds
    % How close must subjects be to marker to start trial?
    marker_dist_threshold = 1; %cm
    % How close must subjects be to screen to count as a touch?
    screen_y_dist_threshold = .5; %cm
    % Space beyond the shape itself that subject can touch to count as a response
    extraSpace = 30;
    % How close must subjects be to target, on both x and y dimensions, to count as a response?
    target_x_dist_threshold = extraSpace*1.5;
    target_y_dist_threshold = extraSpace*2.5;
    
    %% experiment settings
    
    nTrialPerCnd = 24/2; % 24 total over 2 blocks, 12 per block
    nFlanker = 3; % +/-/0
    nBaseOri = 2; % 0/90
    nTarg = 2; % same/diff
    nCrowding = 2; % actual num of crowding conditions; 12 cnd total but 8 who get 12 each
    nCnd = nFlanker*nTarg*nBaseOri;
    
    trialNumber = nTrialPerCnd*nBaseOri*nTarg*nCrowding;
    
    percCatch = 1/3; % 25% catch trials with display in left hemifield (33% of non-catch trials)
    nCatch = trialNumber*percCatch;
    nTotalTrials = nCatch + trialNumber;
    
    if ~resume
        mixedTrialIdx = randperm(nTotalTrials);
        catchIdx = mixedTrialIdx(1:nCatch);
        nonCatchIdx = mixedTrialIdx(nCatch+1:end);

        catchCnd = [repmat(1:4,1,nTrialPerCnd*percCatch/2),repmat(5:8,1,nTrialPerCnd*percCatch),...
            repmat(9:12,1,nTrialPerCnd*percCatch/2)];
        x = randperm(nCatch);
        catchCnd = catchCnd(x);

        nonCatchCnd = [repmat(1:4,1,nTrialPerCnd/2),repmat(5:8,1,nTrialPerCnd),repmat(9:12,1,nTrialPerCnd/2)];
        x = randperm(nTotalTrials - nCatch);
        nonCatchCnd = nonCatchCnd(x);

        cndList = zeros(1,nTotalTrials);
        cndList(catchIdx) = catchCnd;
        cndList(nonCatchIdx) = nonCatchCnd;

        targList = [0 1]; % 0 = same, 1 = diff
        baseOriList = [0 90];
        flankerList = [-1 0 1];
        cnds = CombVec(targList,baseOriList,flankerList);
    end
    
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
    RTdeadline = 10;
    extraMeasurementTime = .2; 
  
    if ~resume
        hemiIndex = ones(1, nTotalTrials);
        hemiIndex(catchIdx) = -1;

        nDiffTrials = sum(mod(cndList,2));
        targTiltList = [-1 1];
        targTiltIdx = repmat(targTiltList,1,nDiffTrials/2);
        n = randperm(nDiffTrials);
        targTiltIdx = targTiltIdx(n);
    end
    
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
    if ~resume
        trials = 0; % initiate trial counter
        difftrials = 0;
        initTime(1:nTotalTrials) = zeros;
        S1_onset_time(1:nTotalTrials) = zeros;
        ISI_onset_time(1:nTotalTrials) = zeros;
        S2_onset_time(1:nTotalTrials) = zeros;
        Mask_onset_time(1:nTotalTrials) = zeros;
        targAngle(1:nTotalTrials) = zeros;
        acc(1:nTotalTrials) = zeros;
        acc_grasp(1:nTotalTrials) = zeros;
        rt(1:nTotalTrials) = zeros;
        Keyresponse(1:nTotalTrials) = zeros;
        rspKey = zeros(1,nTotalTrials);
        rspRatio = [0 0]; % rspRatio(1) counts # lefts, rspRatio(2) counts # rights
    end
    
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
        
        % Subject reached somwhere?
        reachedTo = 0;
        % Stop looking for a reach movement?
        exit = 0;
        
        fprintf(sprintf('OLD TRIAL NUMBER %d',trials))
        trials = trials + 1 
        fprintf(sprintf('NEW TRIAL NUMBER %d',trials))
        
        curCnd = cndList(trials); 
        targID(trials) = cnds(1,curCnd);
        baseOri(trials) = cnds(2,curCnd);
        flanker(trials) = cnds(3,curCnd);
        
        clear pos 
        pos(1,:) = [xCen + eccPx*hemiIndex(trials), yCen]; % center bar
        pos(2,:) = [xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2)]; % UR
        pos(3,:) = [xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2)]; % LL
        pos(4,:) = [xCen + eccPx*hemiIndex(trials) - tfDistPx/sqrt(2), yCen - tfDistPx/sqrt(2)]; % UL
        pos(5,:) = [xCen + eccPx*hemiIndex(trials) + tfDistPx/sqrt(2), yCen + tfDistPx/sqrt(2)]; % LR
        
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
        
        if grasping
            % Reset data holders - these store reach data for each trial
            SOT_data = [];
            xy1_data1 = [];
            xy1_data2 = [];
            currXYZ1_data1 = [];
            currXYZ1_data2 = [];
            currXYZ1_data3 = [];
            xy2_data1 = [];
            xy2_data2 = [];
            currXYZ2_data1 = [];
            currXYZ2_data2 = [];
            currXYZ2_data3 = [];
            currFrame_data = [];

            % Read the tracker once
            clear data
            data = tracker2(5,160);

        %    Wait until they put their finger at the marker to start trial    
            while (abs(marker_pos1(1) - data(3,1)) > marker_dist_threshold ||...
                abs(marker_pos1(2) - data(4,1)) > marker_dist_threshold ||...
                abs(marker_pos1(3) - data(5,1)) > marker_dist_threshold)
                clear data;
                data = tracker2(5,160);
            end
        end
        
        % show fixation cross
        Screen('FillRect', w, grey);
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        initTime(trials) = Screen('Flip',w);
        fixTime = tic;
        if grasping
            % Loop to ensure they stay at the marker for however long "markerWait" is
            while (toc(fixTime) < markerWait(trials)-secPerFrame)
                %clear curr_xyz curr_frame;
                clear data;
                %[curr_xyz1,curr_frame1] = tracker(5,160,1);   
                data = tracker2(5,160); 
                    if (abs(marker_pos1(1) - data(3,1)) > marker_dist_threshold ||...
                        abs(marker_pos1(2) - data(4,1)) > marker_dist_threshold ||...
                        abs(marker_pos1(3) - data(5,1)) > marker_dist_threshold)
                        markerTime = tic;
                    end       
            end 
            old_frame=data(2,1);
        else
            WaitSecs(markerWait(trials));
        end
        
        % draw stimulus display
        for i = 1:5 % n bars
            dstRects(:,i) = CenterRectOnPoint(texrect, pos(i,1), pos(i,2));
        end
        
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
        S1_onset_time(trials) = Screen('Flip',w); % Mark the time when the display went up
        timeElapsed = 0;
        if grasping
            S1_onset = tic;
            while (toc(S1_onset) < T1-secPerFrame)
                %'IN WHILE LOOP'
                % Read the tracker
                clear data;
                %clear curr_xyz curr_frame;
                %[curr_xyz1,curr_frame1] = tracker(5,160,1);
                data = tracker2(5,160);

                % If it's a new sample...
                if data(2,1)~=old_frame
                    %re-define the old frame as current frame
                    old_frame=data(2,1); 
                    %'IN IF STATEMENT!'
                    %Store reach data in the variables
                    % get x and y position of hand in pixel space, and get distance
                    % to screen for real-time feedback
                    xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    curr_pixel_xy1 = [xy1(:,1)];
                    curr_pixel_xy2 = [xy2(:,1)];
                    screen_y_dist1 = norm(data(4,1)-screen_y1);
                    %screen_y_dist2 = norm(data(4,2)-screen_y2);
                    % how much time since stimulus onset?
                    SOT_data = [SOT_data;toc(S1_onset)];
                    % x and y hand positions in pixel space
                    xy1_data1 = [xy1_data1;xy1(1)];
                    xy1_data2 = [xy1_data2;xy1(2)];
                    xy2_data1 = [xy2_data1;xy2(1)];
                    xy2_data2 = [xy2_data2;xy2(2)];
                    % x,y, and z hand positions in real space (cm) from big magnet
                    currXYZ1_data1 = [currXYZ1_data1;data(3,1)];
                    currXYZ1_data2 = [currXYZ1_data2;data(5,1)];
                    currXYZ1_data3 = [currXYZ1_data3;data(4,1)];
                    currXYZ2_data1 = [currXYZ2_data1;data(3,2)];
                    currXYZ2_data2 = [currXYZ2_data2;data(5,2)];
                    currXYZ2_data3 = [currXYZ2_data3;data(4,2)];
                    % What is the current "frame" number (frame meaning number in
                    % sequence of samples recorded by tracker)
                    currFrame_data = [currFrame_data;data(2,1)];
                end
            end
        else
            WaitSecs(T1);
        end
        
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        ISI_onset_time(trials) = Screen('Flip',w,S1_onset_time(trials)+T1-secPerFrame);
        
        if grasping
            ISI_onset = tic;
            while (toc(ISI_onset) < ISI-secPerFrame)
                %'IN WHILE LOOP'
                % Read the tracker
                clear data;
                %clear curr_xyz curr_frame;
                %[curr_xyz1,curr_frame1] = tracker(5,160,1);
                data = tracker2(5,160);

                % If it's a new sample...
                if data(2,1)~=old_frame
                    %re-define the old frame as current frame
                    old_frame=data(2,1); 
                    %'IN IF STATEMENT!'
                    %Store reach data in the variables
                    % get x and y position of hand in pixel space, and get distance
                    % to screen for real-time feedback
                    xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    curr_pixel_xy1 = [xy1(:,1)];
                    curr_pixel_xy2 = [xy2(:,1)];
                    screen_y_dist1 = norm(data(4,1)-screen_y1);
                    %screen_y_dist2 = norm(data(4,2)-screen_y2);
                    % how much time since stimulus onset?
                    SOT_data = [SOT_data;toc(S1_onset)];
                    % x and y hand positions in pixel space
                    xy1_data1 = [xy1_data1;xy1(1)];
                    xy1_data2 = [xy1_data2;xy1(2)];
                    xy2_data1 = [xy2_data1;xy2(1)];
                    xy2_data2 = [xy2_data2;xy2(2)];
                    % x,y, and z hand positions in real space (cm) from big magnet
                    currXYZ1_data1 = [currXYZ1_data1;data(3,1)];
                    currXYZ1_data2 = [currXYZ1_data2;data(5,1)];
                    currXYZ1_data3 = [currXYZ1_data3;data(4,1)];
                    currXYZ2_data1 = [currXYZ2_data1;data(3,2)];
                    currXYZ2_data2 = [currXYZ2_data2;data(5,2)];
                    currXYZ2_data3 = [currXYZ2_data3;data(4,2)];
                    % What is the current "frame" number (frame meaning number in
                    % sequence of samples recorded by tracker)
                    currFrame_data = [currFrame_data;data(2,1)];
                end
            end    
        else
            WaitSecs(ISI);
        end
        
        Screen('DrawTextures', w, barTexVert, [], dstRects, rotAngles2);
        Screen('FillOval', w, stimColor,FIXATION_POSITION,10);
        S2_onset_time(trials) = Screen('Flip',w,ISI_onset_time(trials)+ISI-secPerFrame);
       
        if grasping
            S2_onset = tic;
            while (toc(S2_onset) < T2-secPerFrame)
                %'IN WHILE LOOP'
                % Read the tracker
                clear data;
                %clear curr_xyz curr_frame;
                %[curr_xyz1,curr_frame1] = tracker(5,160,1);
                data = tracker2(5,160);

                % If it's a new sample...
                if data(2,1)~=old_frame
                    %re-define the old frame as current frame
                    old_frame=data(2,1);
                    %'IN IF STATEMENT!'
                    %Store reach data in the variables
                    % get x and y position of hand in pixel space, and get distance
                    % to screen for real-time feedback
                    xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    curr_pixel_xy1 = [xy1(:,1)];
                    curr_pixel_xy2 = [xy2(:,1)];
                    screen_y_dist1 = norm(data(4,1)-screen_y1);
                    %screen_y_dist2 = norm(data(4,2)-screen_y2);
                    % how much time since stimulus onset?
                    SOT_data = [SOT_data;toc(S1_onset)];
                    % x and y hand positions in pixel space
                    xy1_data1 = [xy1_data1;xy1(1)];
                    xy1_data2 = [xy1_data2;xy1(2)];
                    xy2_data1 = [xy2_data1;xy2(1)];
                    xy2_data2 = [xy2_data2;xy2(2)];
                    % x,y, and z hand positions in real space (cm) from big magnet
                    currXYZ1_data1 = [currXYZ1_data1;data(3,1)];
                    currXYZ1_data2 = [currXYZ1_data2;data(5,1)];
                    currXYZ1_data3 = [currXYZ1_data3;data(4,1)];
                    currXYZ2_data1 = [currXYZ2_data1;data(3,2)];
                    currXYZ2_data2 = [currXYZ2_data2;data(5,2)];
                    currXYZ2_data3 = [currXYZ2_data3;data(4,2)];
                    % What is the current "frame" number (frame meaning number in
                    % sequence of samples recorded by tracker)
                    currFrame_data = [currFrame_data;data(2,1)];
                end
            end
        else
            WaitSecs(T2);
        end

        % NOISE CIRCLE
        circRect = SetRect(0,0, barLenPx, barLenPx);       
        for i = 1:5 % n bars
            cDstRects(:,i) = CenterRectOnPoint(circRect, pos(i,1), pos(i,2));
        end       
        if ~flanker(trials)
            cDstRects = cDstRects(:,1);
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
        Mask_onset_time(trials) = Screen('Flip',w,S2_onset_time(trials)+T2-secPerFrame);
        
        if grasping
            if baseOri(trials) % IF HORIZONTAL
                pos_1 = [pos(1,1)+barLenPx/2, pos(1,2)];
                pos_2 = [pos(1,1)-barLenPx/2, pos(1,2)];
            else % IF VERTICAL
                pos_1 = [pos(1,1), pos(1,2)-barLenPx/2];
                pos_2 = [pos(1,1), pos(1,2)+barLenPx/2];
            end
            mask_onset = tic;
            while (toc(mask_onset) < RTdeadline && ~exit)
                %'IN WHILE LOOP'
                % Read the tracker
                clear data;
                data = tracker2(5,160);

                % If it's a new sample...
                if data(2,1)~=old_frame
                    %re-define the old frame as current frame
                    old_frame=data(2,1);
                    %'IN IF STATEMENT!'
                    %Store reach data in the variables
                    % get x and y position of hand in pixel space, and get distance
                    % to screen for real-time feedback
                    xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    curr_pixel_xy1 = [xy1(:,1)];
                    curr_pixel_xy2 = [xy2(:,1)];
                    screen_y_dist1(trials) = norm(data(4,1)-screen_y1);
                    %screen_y_dist2 = norm(data(4,2)-screen_y2);
                    % how much time since stimulus onset?
                    SOT_data = [SOT_data;toc(S1_onset)];
                    % x and y hand positions in pixel space
                    xy1_data1 = [xy1_data1;xy1(1)];
                    xy1_data2 = [xy1_data2;xy1(2)];
                    xy2_data1 = [xy2_data1;xy2(1)];
                    xy2_data2 = [xy2_data2;xy2(2)];
                    % x,y, and z hand positions in real space (cm) from big magnet
                    currXYZ1_data1 = [currXYZ1_data1;data(3,1)];
                    currXYZ1_data2 = [currXYZ1_data2;data(5,1)];
                    currXYZ1_data3 = [currXYZ1_data3;data(4,1)];
                    currXYZ2_data1 = [currXYZ2_data1;data(3,2)];
                    currXYZ2_data2 = [currXYZ2_data2;data(5,2)];
                    currXYZ2_data3 = [currXYZ2_data3;data(4,2)];
                    % What is the current "frame" number (frame meaning number in
                    % sequence of samples recorded by tracker)
                    currFrame_data = [currFrame_data;data(2,1)];

                    % If observers got close enough to the screen, and if they got
                    % within range of a stimulus, end trial and say where they
                    % reached, record time elapsed, and set exit to 1 to exit the
                    % hand recording loop
                    %'Screen touched?'
                    cc = screen_y_dist1(trials) < screen_y_dist_threshold;
                    if screen_y_dist1(trials) < screen_y_dist_threshold && ~exit
                        fprintf(sprintf('Trial %d: Screen touched!\n',trials))
                        %'Screen touched!'
                        reachedTo = 4;
                        if (abs(xy1(1)-pos_1(1))<target_x_dist_threshold && abs(xy1(2)-pos_1(2))...
                                <target_y_dist_threshold) && (abs(xy2(1)-pos_2(1))<target_x_dist_threshold ...
                                && abs(xy2(2)-pos_2(2))<target_y_dist_threshold)
                            reachedTo = 1;
                            fprintf(sprintf('Trial %d: Screen touched & correct grasp!\n',trials))
%                             Beeper(1500,.9,.5); % super high beep
                            acc_grasp(trials) = 1; 
                            timeElapsed = toc(S1_onset);
                        else
                            fprintf(sprintf('Trial %d: Screen touched & incorrect grasp!\n',trials))
                        end
                        if timeElapsed > RTdeadline
                            Beeper(500,.9); % two low beeps if you didn't make the time limit
                            WaitSecs(.3);
                            Beeper(500,.9);
                            exit=1;
                        else
                            Beeper(1000,.9); % high beep if you made the time limit
                            exit=1;
                        end
                    else
                    end

                end
            end
            if ~reachedTo
                Beeper(500,.9);
                WaitSecs(.3);
                Beeper(500,.9);
                timeElapsed = 0;
                reachTime = 0;
                acc_grasp(trials) = 0;
            end
             % mark the time when the end loop begins
            trialLoopEnd = tic;
            % Record movement for an extra "extraMeasurementTime" seconds to get velocity profiles
            while (toc(trialLoopEnd) < extraMeasurementTime)
                % Query the tracker
                clear data;
                data = tracker2(5,160);
                curr_frame = data(2,1);
                if curr_frame~=old_frame
                    %re-define the old frame as current frame
                    old_frame=curr_frame;
                    % Write data to file
                    xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    SOT_data = [SOT_data;toc(S1_onset)];
                    xy1_data1 = [xy1_data1;xy1(1)];
                    xy1_data2 = [xy1_data2;xy1(2)];
                    xy2_data1 = [xy2_data1;xy2(1)];
                    xy2_data2 = [xy2_data2;xy2(2)];
                    currXYZ1_data1 = [currXYZ1_data1;data(3,1)];
                    currXYZ1_data2 = [currXYZ1_data2;data(5,1)];
                    currXYZ1_data3 = [currXYZ1_data3;data(4,1)];
                    currXYZ2_data1 = [currXYZ2_data1;data(3,2)];
                    currXYZ2_data2 = [currXYZ2_data2;data(5,2)];
                    currXYZ2_data3 = [currXYZ2_data3;data(4,2)];
                    currFrame_data = [currFrame_data;data(2,1)];
                    grasp_width = data(3,1) - data(3,2);
                end
            end
        end
       
        keypressed = 0;
        while ~keypressed
            [keyIsDown, secs, responseKey] = KbCheck;
            if keyIsDown
                if responseKey(KbName('q'))
                    ReadPnoRTAllML_ver4(0);
                    Screen('CloseAll');
                    ShowCursor;
                end
                if responseKey(KbName('1'))
                    rspKey(trials) = 0;
                elseif responseKey(KbName('2'))
                    rspKey(trials) = 1;
                end
                if (~targID(trials) && responseKey(KbName('1'))) || (targID(trials) && responseKey(KbName('2'))) 
                    % 1 for same, 2 for diff
                    acc(trials) = 1;
                    Beeper(1000);
                else
                    acc(trials) = 0;
                    Beeper(500);
                end
                rt(trials) = (secs-S1_onset_time(trials));
                Keyresponse(trials) = find(responseKey);
                KbReleaseWait;
                keypressed = 1;
            end
        end
        Screen('FillRect', w, grey);
        Screen('Flip',w);
        
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
        
        if grasping
            save(filenameMatAll);
            cd(oripath);
        end
    end
    
    % Inform subjects that experiment is over
    endDisplay = ['The block is over.\n\n Please go inform the experimenter.'];
    DrawFormattedText(w, endDisplay, 'center', 'center', [255 255 255]);
    Screen('Flip', w);
    WaitSecs(5);

    Screen ('CloseAll');
    ReadPnoRTAllML_ver4(0);
    ShowCursor;
    
    %--------- DONE --------------------
    % rudimentary data analysis
    if blockNum
        ListenChar(1); % Turn keyboard output to command window on
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
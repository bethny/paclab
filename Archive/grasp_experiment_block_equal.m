
%% subject information %%%%%%%%%%%%%%%%%%%
%Screen('Preference', 'SkipSyncTests', 1)
clear all;
KbName('UnifyKeyNames')

VIEWING_DISTANCE_CM = 52;%57;
MONITOR_WIDTH_CM = 44;%35;

WhichScreen = max(Screen('Screens'));
white = WhiteIndex(WhichScreen);
grey = GrayIndex(WhichScreen);
black = BlackIndex(WhichScreen);

subjNum = input('\n Enter subject number: ');

%Create a directory for this subject's data
oripath=pwd;
pathdata=strcat(pwd,filesep,'Subject_folders',filesep,'grasping_',subjNum,filesep);
mkdir(pathdata);
cd(pathdata);
pathdata = pwd;

currBlock = input('\n Enter block number (0=Practice, 1=Exp): ');
resume = input('\n Resume a crashed experiment? (0=no, 1=yes): ');
newBlock = input('\n Start the next block? (0=no, 1=yes): ');
sth=input('\n Enter threshold: ');

calibrate = input('\n Calibrate? ');
calibFile = strcat(pathdata,filesep,subjNum,'_calibration.mat');  %% calibration filename
% if calibrate == 1

cd(oripath);
if exist(calibFile, 'file') ~= 2 || calibrate == 1
    [affine_alignment1,screen_y1,marker_pos1] = calibrate_sensor(1);
    [affine_alignment2,screen_y2,marker_pos2] = calibrate_sensor(2);
    save(calibFile,'affine_alignment1','screen_y1','marker_pos1','affine_alignment2','screen_y2','marker_pos2');
else
    load(calibFile);
end


cd(oripath);% affine_alignment1 = affine_alignment;
% screen_y1 = screen_y;
% marker_pos1 = marker_pos;

%% open window %%%%%%%%%%%%1%%%%%%%%%%%%%%%%%%%

[w, winRect] = Screen('OpenWindow',WhichScreen,128);
MyCLUT = load('gammaTable1.mat');
Screen('LoadNormalizedGammaTable', w, MyCLUT.gammaTable1*[1 1 1]);
[xCen, yCen] = RectCenter(winRect);
[swidth, sheight]=Screen('WindowSize', WhichScreen);
screenCenter = [xCen yCen];
screenInfo.center=screenCenter;
degPerPixel = atan((MONITOR_WIDTH_CM/2)/VIEWING_DISTANCE_CM) * (180/pi) * (2 / swidth);
pixelPerDeg = 1 / degPerPixel;
screenInfo.PPD = pixelPerDeg;
pixelPerCM = swidth / MONITOR_WIDTH_CM;
HideCursor;
secPerFrame = Screen('GetFlipInterval',w);
%% initial value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fixationsize = 0.25;  %size of fixation
size=10;      %size of the grating
vsize=4;      %size of the visible grating
reachsize=6;  %size of the reachable grating
degree=8.5;    %distance between center and grating
key1=97;     %key for 1
key2=98;     %key for 2
sf=1;       %spatial frequency
contrast=1;   %contrast
orientation1=[-1.5*sth 0 1.5*sth;
              -2*sth 0 2*sth];
s=2;
orin1=3;
orientation2=[45 135];
orin2=2;     %the number of orientation2
trialNumber=32;
time1=0.1;  %duration of the 1st grating
time2=0.2;   %duration of ISI
responselimit=1.5;
markerWaitList = [1.75, 2, 2.25];
mn=3;     %the number of markerWait;
gocueList={'grasping','pointing','no action'};
gn=3;     %the number of gocue
pn=2;     %the number of position
rorder=[1,2,2,3,4,5,5,6,7,8,8,9,10,11,11,12];
ron=16;

fdistance=fixationsize*pixelPerDeg;
FIXATION_POSITION = [xCen - fdistance , yCen - fdistance, xCen + fdistance, yCen + fdistance];
distance=degree*pixelPerDeg ; 
x_d=distance;
y_d=distance;
r=reachsize./2*pixelPerDeg;

gaborDimPix = round(size*pixelPerDeg);
VgaborDimPix = round(vsize*pixelPerDeg);
sigma = gaborDimPix/8;
numCycles = size*sf;
freq = numCycles / gaborDimPix;     %spatial frequency of the grating

if currBlock == 0
    nrun=3;
    trialNumber=16;
else
    nrun=12;
end

gocue=[];
for k=1:nrun/3
    gocueIndex=repmat(gocueList,1,3/gn);
    n=randperm(3);
    gocue=[gocue,gocueIndex(1,n)];
end
markerWaitIndex = repmat(markerWaitList,nrun,floor((trialNumber+mn)./mn));
for q=1:nrun
   rn=randperm(trialNumber);
   markerWait(q,:) = markerWaitIndex(q,rn);
end;
medGrey = [150 150 150];
darkGrey = [135 135 135];


% pos(1,:)=[swidth/2-x_d, sheight/2-y_d];
% pos(2,:)=[swidth/2+x_d, sheight/2-y_d];
pos(1,:)=[swidth/2-x_d, sheight/2];
pos(2,:)=[swidth/2+x_d, sheight/2];
pos(3,:)=[swidth/2+x_d, sheight/2+y_d];
pos(4,:)=[swidth/2-x_d, sheight/2+y_d];

orderList=1:ron;
orderIndex=repmat(orderList,nrun,trialNumber./ron);
for q=1:nrun
   n=randperm(trialNumber);
   ordIndex(q,:) = orderIndex(q,n);
end;
         
dstIndex(:,:)=rorder(ordIndex(:,:))>6; 
dstIndex=dstIndex+1;
for q=1:nrun
    for qq=1:trialNumber
        if (rorder(ordIndex(q,qq))==1)||(rorder(ordIndex(q,qq))==4)||(rorder(ordIndex(q,qq))==7)||...
                (rorder(ordIndex(q,qq))==10)
            oriIndex1(q,qq)=1;
        elseif (rorder(ordIndex(q,qq))==2)||(rorder(ordIndex(q,qq))==5)||(rorder(ordIndex(q,qq))==8)||...
                (rorder(ordIndex(q,qq))==11)
            oriIndex1(q,qq)=2;
        else
            oriIndex1(q,qq)=3;
        end;
        if (rorder(ordIndex(q,qq))==1)||(rorder(ordIndex(q,qq))==2)||(rorder(ordIndex(q,qq))==3)||...
                (rorder(ordIndex(q,qq))==7)||(rorder(ordIndex(q,qq))==8)||(rorder(ordIndex(q,qq))==9)
            oriIndex2(q,qq)=1;
        else
            oriIndex2(q,qq)=2;
        end
    end
end

gocue=[];
for k=1:nrun/3
    gocueIndex=repmat(gocueList,1,3/gn);
    n=randperm(3);
    gocue=[gocue,gocueIndex(1,n)];
end

if currBlock==1
    sizeIndex(1:3,1:4)=0;
    for qs=1:3
        sizeIndex2=[];
        sizeList=1:s;
        for k=1:nrun/6
            sizeIndex1 = repmat(sizeList,3,nrun/12);
            n=randperm(2);
            sizeIndex2=[sizeIndex2,sizeIndex1(1,n)];
        end
        sizeIndex(qs,:)=sizeIndex2;
    end
else
    orientation1=[-1.5*sth 0 1.5*sth];
    sizeIndex(1:3,1)=[1; 1; 1];
    gocue={'pointing','grasping','no action'};
end


%%%%%%%%%%%%%%%%%%%%%%
% Tracker Thresholds %
%%%%%%%%%%%%%%%%%%%%%%

% How close must subjects be to marker to start tria?
marker_dist_threshold = 1; %cm
% How close must subjects be to screen to count as a touch?
screen_y_dist_threshold = .2; %cm
% Space beyond the shape itself that subject can touch to count as a
% response
extraSpace = 30;
%shapeRadius = 75;
% How close must subjects be to target, on both x and y dimensions, to
% count as a response?
target_x_dist_threshold = extraSpace* 1.5;
target_y_dist_threshold = extraSpace * 2.5;

chime = MakeBeep(600,.2);
wrongChime = MakeBeep(300,.2);

%%%%%%%%%%%%%%%%%%%%%%
%   Timing variables %
%%%%%%%%%%%%%%%%%%%%%%

% set RTDeadline i.e. how long do they have to make a response
RTDeadline = 10;
% How long to wait in between trials
iti = 1;
% How long to wait for subjects to be at marker before presenting stimuli
%markerWait = 1;
%markerWait = shuffle([ones(numTrials/3,1)*.4;ones(numTrials/3,1)*.5;ones(numTrials/3,1)*.6]);
extraMeasurementTime = .2;

% Preallocate accuracy variable
acc = zeros(nrun,trialNumber);
keyacc = zeros(nrun,trialNumber);
KeyresponseTime=zeros(nrun,trialNumber);
Keyresponse=zeros(nrun,trialNumber);
c(1:3)=0;
task(1:nrun)=0;


if resume == 1 
    load(strcat(pathdata,'\',subjNum,'_',num2str(currBlock),'_grasp_MATDATA'));
    rangeBegin1 = bb;
    rangeBegin2 = i;
    c(task(bb))=c(task(bb))-1;
elseif resume==0 && newBlock==1
    load(strcat(pathdata,'\',subjNum,'_',num2str(currBlock),'_grasp_MATDATA'));
    rangeBegin1 = bb+1;
    rangeBegin2 = 1;
elseif resume==0 && newBlock==0    
    rangeBegin1 = 1;
    rangeBegin2 = 1;
end

m1 = grey + grey*contrast* makegabor(gaborDimPix,freq,sigma,45,VgaborDimPix);
m2 = grey + grey*contrast* makegabor(gaborDimPix,freq,sigma,135,VgaborDimPix);

grasp_Instruction_1('NaN',currBlock, w,m1,m2);

%% main program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for bb=rangeBegin1:nrun
    srcRect=[];
    %Display Instructions
    Screen(w, 'FillRect', grey);
    Screen('Flip', w); 
    grasp_Instruction_1(gocue(bb),currBlock, w,m1,m2);
    
    Screen('TextSize',w, 24);
    if strcmp(gocue(bb),'grasping')
        
        task(bb)=1;
    elseif strcmp(gocue(bb),'pointing')
        
        task(bb)=2;
    elseif strcmp(gocue(bb),'no action')
        
        task(bb)=3;
    end
    c(task(bb))=c(task(bb))+1;
        
    Screen('Flip',w);
    KbWait;
    while KbCheck; end;

    if bb>rangeBegin1
        rangeBegin2=1;
    end
    for i = rangeBegin2:trialNumber

        % For each trial, set the following paramters to zero
        % Subject reached somwhere?
        reachedTo = 0;
        % Stop looking for a reach movement?
        exit = 0;

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
        %clear curr_xyz curr_frame;
        clear data
        data = tracker2(5,160);
        %[curr_xyz2,curr_frame2] = tracker(5,160,2);

    %    Wait until they put their finger at the marker to start trial    
        while (abs(marker_pos1(1) - data(3,1)) > marker_dist_threshold ||...
            abs(marker_pos1(2) - data(4,1)) > marker_dist_threshold ||...
            abs(marker_pos1(3) - data(5,1)) > marker_dist_threshold)
            %clear curr_xyz1 curr_frame1;
            clear data;
            %[curr_xyz1,curr_frame1] = tracker2(5,160,1);
            data = tracker2(5,160);
        end
        Screen('FillRect', w, grey);
        % Draw fixation to indicate the start of the trial       
        Screen('FillOval', w, white,FIXATION_POSITION,10);     
        FixTime(bb,i)=Screen('Flip',w);

        % Loop to ensure they stay at the marker for however long "markerWait" is
        markerTime = tic;
        while (toc(markerTime) < markerWait(bb,i)-secPerFrame)
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
        %curr_frame = data(2,1);
        old_frame=data(2,1);
        
        
        r1=orientation1(sizeIndex(task(bb),c(task(bb))),oriIndex1(bb,i));
        r2=orientation2(oriIndex2(bb,i));
        if sizeIndex(task(bb),c(task(bb)))==1
            ss(bb)=1;
        else
            ss(bb)=2;
        end
        %make Gabor 1
        m = grey + grey*contrast* makegabor(gaborDimPix,freq,sigma,r2+r1,VgaborDimPix);
        gabortex=Screen('MakeTexture', w, m);
        texrect = Screen('Rect', gabortex);
        dstRects = CenterRectOnPoint(texrect,pos(dstIndex(bb,i),1), pos(dstIndex(bb,i),2))';        
        Screen('DrawTexture', w, gabortex, srcRect,dstRects);    
        Screen('FillOval', w,white,FIXATION_POSITION,10);  
        
        %make Gabor 2
        m2 = grey + grey*contrast* makegabor(gaborDimPix,freq,sigma,r2,VgaborDimPix);
        gabortex2=Screen('MakeTexture', w, m2);
        texrect2 = Screen('Rect', gabortex2);
        dstRects2 = CenterRectOnPoint(texrect2,pos(dstIndex(bb,i),1), pos(dstIndex(bb,i),2))';
        FirstShow(bb,i)=Screen('Flip',w);
        
        FirstShowtic(bb,i)=tic;
        timeElapsed=0;
        while (toc(FirstShowtic(bb,i)) < time1-secPerFrame)
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
                SOT_data = [SOT_data;toc(FirstShowtic(bb,i))];
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
            end;
        end;
        Screen('FillOval', w, white,FIXATION_POSITION,10);          
        ISITime(bb,i)=Screen('Flip',w,FirstShow(bb,i)+time1-secPerFrame);
        ISITimetic(bb,i)=tic;
        while (toc(ISITimetic(bb,i)) < time2-secPerFrame)
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
                SOT_data = [SOT_data;toc(FirstShowtic(bb,i))];
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
            end;
        end;        
        Screen('DrawTexture', w, gabortex2, srcRect,dstRects2);
        %Screen('PutImage', w, m,dstRects);
        SecondShow(bb,i)=Screen('Flip',w,ISITime(bb,i)+time2-secPerFrame);
        
        % Mark the time when the display went up
        SecondShowtic(bb,i) = tic;
        if oriIndex2(bb,i)==2
            pos_top = [pos(dstIndex(bb,i),1)+r./sqrt(2), pos(dstIndex(bb,i),2)-r./sqrt(2)];
            pos_bottom = [pos(dstIndex(bb,i),1)-r./sqrt(2), pos(dstIndex(bb,i),2)+r./sqrt(2)];
        elseif oriIndex2(bb,i)==1
            pos_top = [pos(dstIndex(bb,i),1)-r./sqrt(2), pos(dstIndex(bb,i),2)-r./sqrt(2)];
            pos_bottom = [pos(dstIndex(bb,i),1)+r./sqrt(2), pos(dstIndex(bb,i),2)+r./sqrt(2)];
        end;

        %'IN WHILE LOOP?'
        % Loop to look for reach movements
        while (toc(SecondShowtic(bb,i)) < RTDeadline && ~exit)
            if strcmp(gocue(bb),'no action')
                exit=1;
                reachedTo = 3;
                acc(bb,i) = 0; 
                timeElapsed = 0;
            end;
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
                SOT_data = [SOT_data;toc(FirstShowtic(bb,i))];
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

                if (abs(marker_pos1(1) - data(3,1)) > marker_dist_threshold ||...
                    abs(marker_pos1(2) - data(4,1)) > marker_dist_threshold ||...
                    abs(marker_pos1(3) - data(5,1)) > marker_dist_threshold)
    %                 reach_onset_time = tic;
                end

                 % If observers got close enough to the screen, and if they got
                % within range of a stimulus, end trial and say where they
                % reached, record time elapsed, and set exit to 1 to exit the
                % hand recording loop
                %'Screen touched?'
                %rect_cen_y
                %rect_cen_y_bottom 
                cc=screen_y_dist1 < screen_y_dist_threshold;
                 if screen_y_dist1 < screen_y_dist_threshold && ~exit
                         %'Screen touched!'
                         reachedTo=4;
                         if strcmp(gocue(bb),'grasping')
                             
                             if (abs(xy1(1)-pos_top(1))<target_x_dist_threshold && abs(xy1(2)-pos_top(2))...
                                     <target_y_dist_threshold) && (abs(xy2(1)-pos_bottom(1))<target_x_dist_threshold ...
                                     && abs(xy2(2)-pos_bottom(2))<target_y_dist_threshold) 
                                reachedTo = 1;
                                acc(bb,i) = 1; 
                                timeElapsed = toc(FirstShowtic(bb,i));
                                
                             end
                         elseif strcmp(gocue(bb),'pointing')
                             if (abs(xy1(1)-pos(dstIndex(bb,i),1))<target_x_dist_threshold && ...
                                     abs(xy1(2)-pos(dstIndex(bb,i),2))<target_y_dist_threshold)
                                reachedTo = 2;
                                acc(bb,i) = 1; 
                                timeElapsed = toc(FirstShowtic(bb,i));
                             end
                         elseif  strcmp(gocue(bb),'no action')
                                reachedTo = 3;
                                acc(bb,i) = 0; 
                                timeElapsed = 0;
                         end
                         if (timeElapsed>responselimit)
                            Snd('Play',wrongChime);
                            WaitSecs(.3);
                            Snd('Play',wrongChime);
                            exit=1;
                         else
                            Snd('Play',chime);
                            exit=1;
                         end
                 end

            end
        end
        if (reachedTo == 0) 
            Snd('Play',wrongChime);
            WaitSecs(.3);
            Snd('Play',wrongChime); 
            timeElapsed = 0;
            reachTime = 0;
            acc(bb,i) = 0;
        end

        % mark the time when the end loop begins
        trialLoopEnd = tic;
        % Record movement for an extra "extraMeasurementTime" seconds to get velocity profiles
            while (toc(trialLoopEnd) < extraMeasurementTime)
                % Query the tracker
                %clear curr_xyz curr_frame;
                clear data;
                %[curr_xyz1,curr_frame1] = tracker(5,160,1);
                data = tracker2(5,160);
                curr_frame = data(2,1);
                if curr_frame~=old_frame
                    %re-define the old frame as current frame
                    old_frame=curr_frame;
                    % Write data to file
                    xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    SOT_data = [SOT_data;toc(FirstShowtic(bb,i))];
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
        keypressed=0;
        while (keypressed==0) 
            [keyIsDown, secs, responseKey]=KbCheck; 

            if keyIsDown
                if responseKey(KbName('ESCAPE'))
                    %Close the program, shutdown tracker
                    ReadPnoRTAllML_ver4(0);
                    ShowCursor;
                    Priority(0);
                    Screen('CloseAll');
                end
                if ((find(responseKey)==key1)&&(r1==0))||((find(responseKey)==key2)&&(r1~=0))
                    keyacc(bb,i)=1;
                else
                    keyacc(bb,i)=0;
                end;
                KbReleaseWait;
                KeyresponseTime(bb,i) = (secs-SecondShowtic(bb,i));
                Keyresponse(bb,i)=find(responseKey);               
                keypressed=1;
                Snd('Play',chime);
            end;

        end

       Screen('FillRect', w, grey);
       Screen('Flip',w);
            % Write this trial's data to two files - exp. data and
            % tracker data
            for x = 1:length(SOT_data)
                dlmwrite(strcat(pathdata,'\',subjNum,'_',num2str(currBlock),'_graspExp.txt'),...
                    [bb,i,SOT_data(x),xy1_data1(x),xy1_data2(x),xy2_data1(x),xy2_data2(x),currXYZ1_data1(x),currXYZ1_data2(x),currXYZ1_data3(x),currXYZ2_data1(x),currXYZ2_data2(x),currXYZ2_data3(x),currFrame_data(x)],'-append', 'roffset', [],'delimiter', '\t');
            end

            %%%FIGURE THIS OUT: OUTPUT SPEED AND WIDTH DATA TO
            %%%EASIER-TO-READ FILE
            dlmwrite(strcat(pathdata,'\',subjNum,'_',num2str(currBlock),'_graspExp_block',num2str(bb),'.txt'),...
                [bb, i, dstIndex(bb,i),r1,r2, reachedTo, ISITime(bb,i)-FirstShow(bb,i), SecondShow(bb,i)-FirstShow(bb,i),timeElapsed, Keyresponse(bb,i),keyacc(bb,i),acc(bb,i),  markerWait(bb,i)],'-append', 'roffset', [],'delimiter', '\t');

            save(strcat(pathdata,'\',subjNum,'_',num2str(currBlock),'_grasp_MATDATA'));  
    end
        % Inform subjects that experiment is over, wait for keypress
    endDisplay = ['The block is over.\n\n Please go inform the experimenter.'];
    DrawFormattedText(w, endDisplay, 'center', 'center', black);
    Screen('Flip', w);
    KbWait;
    while KbCheck; end;
    if currBlock==1
        %Close the program, shutdown tracker
        ReadPnoRTAllML_ver4(0);
        ShowCursor;
        Priority(0);
        Screen('CloseAll');
    end
end
% cd (pathdata);
% cd ..
%Close the program, shutdown tracker
ReadPnoRTAllML_ver4(0);
ShowCursor;
Priority(0);
Screen('CloseAll');




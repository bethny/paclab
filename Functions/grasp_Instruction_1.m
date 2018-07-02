function grasp_Instruction_1(task,currBlock,w,m1,m2)

%==============================================================
% displayInstructions.m
%==============================================================

Screen('FillRect', w, [127 127 127]);
Screen('TextSize', w, []);
[swidth, sheight]=Screen('WindowSize', w);
instText = [];
srcRect1=[];
srcRect2=[];
%Practice Instructions
if strcmp(task,'NaN')
    instText = ['Each trial will begin with a fixation cross.\n\n'...
    'After maintaining fixation of the cross, two targets will appear one by one.\n\n\n\n'...
    'Just after you see the first target,\n\n'...'
    'You should start to grasp/point/perform no action on the second target as soon as possible. \n\n'...
    'If you respond in time, you will hear a beep. \n\n'...
    'If you run out of time, you will hear two beeps. \n\n\n\n'...
    'After that, your need report if the orientations of two targets are same or not. \n\n'...
    'If the orientations are same, push button 1, \n\n'...
    'If the orientations are different, push button 2.'];
    gabortex1=Screen('MakeTexture', w, m1);
    texrect1 = Screen('Rect', gabortex1);
    dstRects1 = CenterRectOnPoint(texrect1,swidth/2-200, sheight/2+300)';  
    Screen('DrawTexture', w, gabortex1, srcRect1,dstRects1);
    gabortex2=Screen('MakeTexture', w, m2);
    texrect2 = Screen('Rect', gabortex2);
    dstRects2 = CenterRectOnPoint(texrect2,swidth/2+200, sheight/2+300)';  
    Screen('DrawTexture', w, gabortex2, srcRect2,dstRects2);    
elseif strcmp(task,'grasping')
    instText = ['Move to grasp the second target after you see the first target.\n\n'... 
                'If the orientations are same, push button 1, \n\n'...
                'If the orientations are different, push button 2. \n\n\n\n'...
                'Please tell the experimenter when you are ready. \n\n'];
elseif strcmp(task,'pointing')
    instText = ['Move to point the second target after you see the first target.\n\n'... 
                'If the orientations are same, push button 1, \n\n'...
                'If the orientations are different, push button 2. \n\n\n\n'...
                'Please tell the experimenter when you are ready. \n\n'];            
elseif strcmp(task,'no action')
    instText = ['Perform no action to the second target.\n\n'... 
                'If the orientations are same, push button 1, \n\n'...
                'If the orientations are different, push button 2. \n\n\n\n'...
                'Please tell the experimenter when you are ready. \n\n'];
end
DrawFormattedText(w, instText, 'Center', 'Center', [0 0 0]);
% if currBlock==5

% end

Screen('Flip', w);
WaitSecs(1);
FlushEvents('keyDown');
GetChar;
Screen('TextSize', w, 25);
Screen('Flip', w);


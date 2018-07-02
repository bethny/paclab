function [affine_alignment,screen_y,marker_pos] = calibrate_sensor(sensorNum)
% function [affine_alignment,screen_y,marker_pos] = calibrate_sensors(sensorNum)
% Calibrate one tracker sensor at a time, eg for grasping studies
% > sensorNum (int): 1 = index finger, 2 = thumb
% < affine_alignment, screen_y, marker_pos = same outputs as calibrate_tracker

[data,~,~] = ReadPnoRTAllML_ver4(5);
    
Snd('Open');
beep = MakeBeep(500,.2);
Screen('Preference','SkipSyncTests',1);
window = Screen('OpenWindow',max(Screen('Screens')));

if sensorNum == 1
    DrawFormattedText(window,'Place index finger on the marker on the table.','center','center');
elseif sensorNum == 2
    DrawFormattedText(window,'Place thumb on the marker on the table.','center','center');
end
Screen(window,'Flip')         

KbWait;
while KbCheck; end;

KbWait;
while KbCheck; end;

[data,~,~] =ReadPnoRTAllML_ver4(5);
marker_pos = data(3:5,sensorNum)';

Screen(window,'FillRect',WhiteIndex(window));
res = Screen('Resolution',window);
rectSize = 20;  
rect = [0 0 rectSize rectSize];
xres = res.width; yres = res.height;
ft_pos = [];
screen_pos = [];
screen_y = [];
numx = 3; numy = 3;

for y = 1:numy
    for x = 1:numx
        target =  CenterRectOnPoint(rect,(x/(numx+1))*xres,(y/(numy+1))*yres);
        Screen(window,'FillRect',BlackIndex(window),target)
        Screen(window,'Flip')     
        screen_pos = [screen_pos,[x/(numx+1)*xres; y/(numy+1)*yres]];          
        
        KbWait;
        while KbCheck; end;

        [data,~,~] =ReadPnoRTAllML_ver4(5);
        pos = data(3:5,sensorNum);
        ft_pos = [ft_pos,pos([1 3])];          
        screen_y = [screen_y pos(2)];
        Snd('Play',beep);
        
    end
end
save('cal.mat');
ReadPnoRTAllML_ver4(0);
DrawFormattedText(window,'Done calibrating this sensor.','center','center');
Screen(window,'Flip')         
WaitSecs(.3)
Snd('Close');

affine_alignment = [];
for i = 1:2
    affine_alignment(i,:) = ([ft_pos' ones(numx*numy,1)]\screen_pos(i,:)');
end
screen_y = sum(screen_y)/length(screen_y); 
save('calibrate');
sca

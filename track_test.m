% try
clear all;

%Screen('Preference','SkipSyncTests',1);
window = Screen('OpenWindow',max(Screen('Screens')));
Screen(window,'FillRect',WhiteIndex(window));
res = Screen('Resolution',window);
Screen('Flip',window,0,1); 
stimrect = [0,0,20,20];
xres = res.width; yres = res.height;
x = xres/2;
y = yres/2;

pathdata=pwd;
h=actxserver('Prok.Liberty');
h.connect;
h.start;
old_frame=0;
SOT_data = [];
currXYZ1_data1 = [];
currXYZ1_data2 = [];
currXYZ1_data3 = [];
currXYZ2_data1 = [];
currXYZ2_data2 = [];
currXYZ2_data3 = [];
currFrame_data = [];
c=1;
s=tic;
for i=1
    h.flush;
    while (toc(s) < 10)
        %t(c)=Screen('Flip',window);
        Screen('FillOval', window ,[1 1 1],CenterRectOnPoint(stimrect,x,y));
        t1(c)=Screen('Flip',window);
        FirstShowtic(c)=tic;
        timeElapsed=0;
        while (toc(FirstShowtic(c)) < 0.1-0.017)
            %[curr_xyz1,curr_frame1] = tracker(5,160,1);
                data = h.position;
                %h.get
                new_frame = h.available_samples;
                % If it's a new sample...
                if new_frame~=old_frame
                    %re-define the old frame as current frame
                    old_frame=new_frame; 

                    %'IN IF STATEMENT!'
                    %Store reach data in the variables
                    % get x and y position of hand in pixel space, and get distance
                    % to screen for real-time feedback
                    %xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    %xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    %curr_pixel_xy1 = [xy1(:,1)];
                    %curr_pixel_xy2 = [xy2(:,1)];
                    %screen_y_dist1 = norm(data(2,1)-screen_y1);
                    %screen_y_dist2 = norm(data(4,2)-screen_y2);
                    % how much time since stimulus onset?
                    SOT_data = [SOT_data;toc(s)];
                    % x,y, and z hand positions in real space (cm) from big magnet
                    currXYZ1_data1 = [currXYZ1_data1;data(1,1)];
                    currXYZ1_data2 = [currXYZ1_data2;data(1,3)];
                    currXYZ1_data3 = [currXYZ1_data3;data(1,2)];
    %                 currXYZ2_data1 = [currXYZ2_data1;data(1,7)];
    %                 currXYZ2_data2 = [currXYZ2_data2;data(1,9)];
    %                 currXYZ2_data3 = [currXYZ2_data3;data(1,8)];
                    % What is the current "frame" number (frame meaning number in
                    % sequence of samples recorded by tracker)
                    currFrame_data = [currFrame_data;new_frame];
                end;
        end;
        
        t2(c)=Screen('Flip',window,t1(c)+0.1-0.017);
        SecondShowtic(c)=tic;
       
        while (toc(SecondShowtic(c)) < 0.2-0.017)
        %disp(toc(s));
                %'IN WHILE LOOP'
                % Read the tracker
                %clear data;
                %clear curr_xyz curr_frame;
                %[curr_xyz1,curr_frame1] = tracker(5,160,1);
                data = h.position;
                %h.get
                new_frame = h.available_samples;
                % If it's a new sample...
                if new_frame~=old_frame
                    %re-define the old frame as current frame
                    old_frame=new_frame; 

                    %'IN IF STATEMENT!'
                    %Store reach data in the variables
                    % get x and y position of hand in pixel space, and get distance
                    % to screen for real-time feedback
                    %xy1 = bsxfun(@plus,affine_alignment1(:,3),affine_alignment1(:,1:2)*[data(3,1);data(5,1)]);
                    %xy2 = bsxfun(@plus,affine_alignment2(:,3),affine_alignment2(:,1:2)*[data(3,2);data(5,2)]);
                    %curr_pixel_xy1 = [xy1(:,1)];
                    %curr_pixel_xy2 = [xy2(:,1)];
                    %screen_y_dist1 = norm(data(2,1)-screen_y1);
                    %screen_y_dist2 = norm(data(4,2)-screen_y2);
                    % how much time since stimulus onset?
                    SOT_data = [SOT_data;toc(s)];
                    % x,y, and z hand positions in real space (cm) from big magnet
                    currXYZ1_data1 = [currXYZ1_data1;data(1,1)];
                    currXYZ1_data2 = [currXYZ1_data2;data(1,3)];
                    currXYZ1_data3 = [currXYZ1_data3;data(1,2)];
    %                 currXYZ2_data1 = [currXYZ2_data1;data(1,7)];
    %                 currXYZ2_data2 = [currXYZ2_data2;data(1,9)];
    %                 currXYZ2_data3 = [currXYZ2_data3;data(1,8)];
                    % What is the current "frame" number (frame meaning number in
                    % sequence of samples recorded by tracker)
                    currFrame_data = [currFrame_data;new_frame];

                     % If observers got close enough to the screen, and if they got
                    % within range of a stimulus, end trial and say where they
                    % reached, record time elapsed, and set exit to 1 to exit the
                    % hand recording loop
                    %'Screen touched?'
                    %rect_cen_y
                end
        end
        Screen('FillOval', window ,[1 0 1],CenterRectOnPoint(stimrect,x,y));
        t3(c)=Screen('Flip',window,t2(c)+0.2-0.017);
        c=c+1;
    end
                
    for x = 1:length(SOT_data)
                dlmwrite(strcat(pathdata,'/graspExp5.txt'),...
                    [currXYZ1_data1(x),currXYZ1_data2(x),currXYZ1_data3(x),],'-append', 'roffset', [],'delimiter', '\t');
                   % [currXYZ1_data1(x),currXYZ1_data2(x),currXYZ1_data3(x),currXYZ2_data1(x),currXYZ2_data2(x),currXYZ2_data3(x)],'-append', 'roffset', [],'delimiter', '\t');
    end
    
end;
Screen('CloseAll');
h.stop;
h.disconnect;

% catch
%     fprintf('Failed' );
%     Screen('CloseAll');
%     h.stop;
%     h.disconnect;
% end
%% STAY DURING PRE-STIM

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

%% STAY DURING S1

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
    end
end

%% STAY DURING ISI

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


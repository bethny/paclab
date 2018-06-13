function point = gridGen(nRow, nCol, upperBound, lowerBound, leftBound, rightBound)

% diff x points
x = linspace(leftBound, rightBound, nCol);

% diff y points
y = linspace(upperBound, lowerBound, nRow);

[gridX, gridY] = meshgrid(x,y); 

idx = 0;
for i = 1:size(gridX,1)
    for j = 1:size(gridX,2)
        idx = idx + 1;
        point(idx,:) = [gridX(i,j), gridY(i,j)];
    end
end

%% DEBUGGING
% colors = [...
%     {hex2rgb('ff0000')},{hex2rgb('FA4E00')},{hex2rgb('F69900')},{hex2rgb('F2E200')},...
%     {hex2rgb('B3ED00')},{hex2rgb('67E900')},{hex2rgb('1DE500')},{hex2rgb('00E128')},...
%     {hex2rgb('00DC6D')},{hex2rgb('00D8AE')},{hex2rgb('00BBD4')},{hex2rgb('0076D0')},...
%     {hex2rgb('0034CB')},{hex2rgb('0A00C7')},{hex2rgb('4700C3')},{hex2rgb('8100BF')}];
% close all
% 
% figure
% hold on
% for i = 1:length(point2)
%     plot(point2(i,1),point2(i,2),'o','color',colors{i})
% end

end
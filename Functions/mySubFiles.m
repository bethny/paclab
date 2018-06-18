function dataFiles = mySubFiles(dir,descriptor,idx)

%   Description:    Lists all the relevant subfiles in a directory, in order. 
%
%   Input:
%           dir:        String. Directory.
%           descriptor: String. Denotes keyword for type of file you want
%                       to find.
%           idx:        Double. Letter index of desired keyword.
%
%   Output:
%           dataFiles:  String array. Lists desired files in numeric order.
    if nargin < 3
        idx = 1;
    end
    allFiles = ls(dir);
    splitFiles = strsplit(allFiles);
    x = 1;
    for i = 1:length(splitFiles)-1
        if cell2mat(strfind(splitFiles(i),descriptor)) == idx
            dataFiles(x) = splitFiles(i);
            x = x+1;
        end
    end
    dataFiles = sort(dataFiles)';

end
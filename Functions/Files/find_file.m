function [path, num, counter] = find_file(str, searchPath)
%% Find file in directory containing specified string or pattern

if nargin < 2
    
    searchPath = pwd;
end

name = [];
num = nan;
counter = 1;

folderData = struct2cell(dir(searchPath));
dirNames = folderData(1,:);

for i = 1:length(dirNames)
    
    if contains(dirNames(i), str)
        
        name = dirNames{i};
        num = i;
        counter = counter + 1;
    end
end

path = fullfile(searchPath, name);
function [path, counter] = find_file(filepattern, search_path)
%FIND_FILE searches for specified file/directory string pattern within 
%directory provided.
%   Inputs:
%   filepattern - String with file/directory name or partial name
%   search_path - Directory to search within 
%       (default = current directoy)
%
%   Outputs:
%   path - path to last file/directory matching filepattern
%   counter - number of matches found, useful for ensuring file/directory
%   names are unique.

if nargin < 2, search_path = pwd; end

name = [];
counter = 1;

% Collecting names
folder_data = struct2cell(dir(search_path));
names = folder_data(1,:);

% Loop through names to find matches
for i = 1:length(names)
    
    if contains(names(i), filepattern)
        
        name = names{i};
        counter = counter + 1;
    end
end

path = fullfile(search_path, name);
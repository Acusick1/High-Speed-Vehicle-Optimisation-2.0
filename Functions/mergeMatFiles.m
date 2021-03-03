function mergeMatFiles(folder, fields)
%MERGEMATFILES loads all mat files in path and concatenates them
%   folder: path to load files from (default = pwd)
%   fields: specific fields to merge, cell array of strings (default = all)

if nargin < 1 || isempty(folder), folder = pwd; end

FileList = dir(fullfile(folder, '*.mat'));  % List of all MAT files

if nargin < 2 || isempty(fields)
    
    temp = load(fullfile(folder, FileList(1).name));
    fields = fieldnames(temp);
end

allData = struct();
for iFile = 1:numel(FileList)              % Loop over found files
    
    data = load(fullfile(folder, FileList(iFile).name), fields{:});
    
    for iField = 1:numel(fields)              % Loop over fields of current file
        
        aField = fields{iField};
        
        if isfield(allData, aField)             % Attach new data:
            
            allData.(aField) = [allData.(aField); data.(aField)];
        else
            allData.(aField) = data.(aField);
        end
    end
end

allData.filenames = {FileList.name}';
save(fullfile(folder, '_AllData.mat'), '-struct', 'allData');
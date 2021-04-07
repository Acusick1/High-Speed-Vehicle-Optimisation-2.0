function mergeMatFiles(folder, inFile, fields, outFile)
%MERGEMATFILES loads all mat files in path and concatenates them
%   folder: path to load files from (default = pwd)
%   fields: specific fields to merge, cell array of strings (default = all)

if nargin < 1 || isempty(folder), folder = pwd; end
if nargin < 2 || isempty(inFile), inFile = '*'; end

FileList = dir(fullfile(folder, [inFile '.mat']));  % List of all MAT files

if nargin < 3 || isempty(fields)
    
    temp = load(fullfile(folder, FileList(1).name));
    fields = fieldnames(temp);
end

allData = struct();
for iFile = 1:numel(FileList)              % Loop over found files
    
    data = load(fullfile(FileList(iFile).folder, FileList(iFile).name), fields{:});
    
    for iField = 1:numel(fields)              % Loop over fields of current file
        
        aField = fields{iField};
        
        if isfield(allData, aField)             % Attach new data:
            
            allData.(aField) = [allData.(aField); data.(aField)];
        else
            allData.(aField) = data.(aField);
        end
    end
end

% Removing wildcards for saving
while contains(folder, '*')
    
    [folder, ~] = fileparts(folder);
end

allData.filenames = {FileList.name}';
save(fullfile(folder, outFile), '-struct', 'allData');
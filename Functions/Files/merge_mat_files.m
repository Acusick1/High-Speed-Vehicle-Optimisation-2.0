function merge_mat_files(folder, file_pattern, fields, output_filename)
%MERGE_MAT_FILES loads mat files and concatenates them.
%   Inputs:
%   folder - path to load files from (default = pwd)
%   file_pattern - pattern to load mat files (default = all files)
%   fields - specific fields to merge, cell array of strings 
%       (default = all)
%   output_filename - filename of merged file

if nargin < 1 || isempty(folder), folder = pwd; end
if nargin < 2 || isempty(file_pattern), file_pattern = '*'; end

% List of all MAT files
file_list = dir(fullfile(folder, [file_pattern '.mat']));

if nargin < 3 || isempty(fields)
    
    temp = load(fullfile(folder, file_list(1).name));
    fields = fieldnames(temp);
end

concat_data = struct();
for i = 1:numel(file_list) % Loop over found files
    
    filename = fullfile(file_list(i).folder, file_list(i).name);
    data = load(filename, fields{:});
    
    % Loop over fields of current file
    for j = 1:numel(fields)
        
        field = fields{j};
        
        % Attach new data
        if isfield(concat_data, field)
            
            concat_data.(field) = [concat_data.(field); data.(field)];
        else
            concat_data.(field) = data.(field);
        end
    end
end

% Removing wildcards for saving
while contains(folder, '*')
    
    [folder, ~] = fileparts(folder);
end

% Saving list of merged files
concat_data.filenames = {file_list.name}';
save(fullfile(folder, output_filename), '-struct', 'concat_data');
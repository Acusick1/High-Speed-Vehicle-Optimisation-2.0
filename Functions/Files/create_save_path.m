function [save_path, save_dir] = create_save_path(save_dir)
%CREATE_SAVE_PATH makes a unique folder based on input within the
%environment's results directory.

results_path = get_results_path();

% Look for existing save_dir, append number to make it unique.
[~, i] = find_file(save_dir, results_path);
save_dir = [save_dir num2str(i)];

save_path = fullfile(results_path, save_dir);

% Appending number is not robust, so make sure folder does not exist
if ~exist(save_path, 'dir')
    mkdir(save_path)
else
    error("Directory '%s' already exists", save_path)
end
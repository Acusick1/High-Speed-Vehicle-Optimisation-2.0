function [save_path, save_dir] = init_save_path(save_dir)

save_path = fullfile(pwd, 'Results');
[~,~,i] = find_file(save_dir, save_path);
save_dir = [save_dir num2str(i)];
save_path = fullfile(save_path, save_dir);
mkdir(save_path)
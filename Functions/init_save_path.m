function save_path = init_save_path(save_dir)

save_path = fullfile(pwd, 'Results');
[~,~,i] = find_file(save_dir, save_path);
save_path = fullfile(save_path, [save_dir num2str(i)]);
mkdir(save_path)
function save_figs(dir)
%SAVE_FIGS saves all open figures to specified directory

fig_list = findobj(allchild(0), 'flat', 'Type', 'figure');

for i = 1:length(fig_list)
  
    fig_handle = fig_list(i);
    fig_name   = get(fig_handle, 'Name');
    
    if isempty(fig_name), fig_name = num2str(i); end
    
    savefig(fig_handle, fullfile(dir, fig_name));
end

close all
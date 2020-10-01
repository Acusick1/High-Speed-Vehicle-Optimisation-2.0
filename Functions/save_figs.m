function save_figs(dir)

figList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(figList)
  
    figHandle = figList(iFig);
    figName   = get(figHandle, 'Name');
    if isempty(figName), figName = num2str(iFig); end
    savefig(figHandle, fullfile(dir, figName));
end

close all
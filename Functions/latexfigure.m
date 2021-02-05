function latexfigure(fig, sym)

if nargin < 1 || isempty(fig), fig = gcf; end
if nargin < 2, sym = [0 0]; end

axs = fig.Children;
xyz = 'xyz';

for i = 1:numel(axs)
    
    ax = axs(i);
    set(ax,'XMinorTick','on','YMinorTick','on','box','on','TickLabelInterpreter','latex')%,...
        %'DataAspectRatio',[1 1 1])
    
    for j = 1:numel(sym)
        
        label = [xyz(j) 'label'];
        str = get(get(gca, label), 'String');
        
        if sym(i), str = ['$' str '$']; end
        
        feval(label, str, 'Interpreter', 'latex');
    end
end
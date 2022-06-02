function latex_figure(fig, symbolic)
%LATEX_FIGURE takes input figure and converts all subplots/axes ticks and
%labels to latex interpreter. 
%   Inputs:
%   fig - Figure handle to convert (default = gcf).
%   symbolic - Boolean array to use latex math mode on symbolic axis label,
%       either 1x2 array for [x, y] or 1x3 [x, y, z] (default = [0 0] for x
%       and y axes).

if nargin < 1 || isempty(fig), fig = gcf; end
if nargin < 2, symbolic = [0 0]; end

axs = fig.Children;
xyz = 'xyz';

for i = 1:numel(axs)
    % Loop through subplots/axes within figure
    ax = axs(i);
    set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on', 'box', 'on', ...
        'TickLabelInterpreter', 'latex')
    
    for j = 1:numel(symbolic)
        % Loop through x, y, and z axis if specified
        label = [xyz(j) 'label'];
        str = get(get(gca, label), 'String');
        
        if symbolic(i)
            feval(label, ['$' str '$'], 'Interpreter', 'latex');
        else
            feval(label, str, 'Interpreter', 'latex');
        end
    end
end
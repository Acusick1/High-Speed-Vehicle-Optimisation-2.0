function plot_format(f, axes_spec, line_spec)
%PLOT_FORMAT formats a given formats to the specifications below
%   Given a figure handle input, each axis, and each line within will be
%   specified according to axesSpec and lineSpec respectively
%   Inputs:
%       f        - figure handle, or current figure if not specified
%       axesSpec - axes formatting, specified as name value pairs within a
%                  cell array, see default value below
%       lineSpec - as above for individual lines

if nargin < 1 || isempty(f), f = gcf; end
if nargin < 2 || isempty(axes_spec)
    
    axes_spec = {'FontSize',            14      , ...
                'box',                  'on'    , ...
                'XMinorTick',           'on'    , ...
                'YMinorTick',           'on'    , ...
                'ZMinorTick',           'on'    , ...
                'TickLabelInterpreter', 'latex'};
end
if nargin < 3 || isempty(axes_spec)
    
    line_spec = {'LineWidth',            1.5    , ...
                'MarkerSize',            18};
end

legend_spec = {'Interpreter',      'latex'};
axs = get(f, 'Children');

% Loop through axes and carry out formatting
for i = 1:numel(axs)
    
    if strcmp(get(axs(i),'tag'), 'legend')
        
        set(axs(i), legend_spec{:});
    
    elseif isempty(get(axs(i),'tag'))
    
        set(axs(i), axes_spec{:})
        
        % Find all data within axes and format if line object
        data = get(axs(i), 'Children');
        for j = 1:numel(data)
            
            if isa(data(j), 'matlab.graphics.chart.primitive.Line')
            
                set(data(j), line_spec{:})
            end
        end
    end
end

end
function plotFormat(f, axesSpec, lineSpec)
%PLOTFORMAT formats a given formats to the specifications below
%   Given a figure handle input, each axis, and each line within will be
%   specified according to axesSpec and lineSpec respectively
%   Inputs:
%       f        - figure handle, or current figure if not specified
%       axesSpec - axes formatting, specified as name value pairs within a
%                  cell array, see default value below
%       lineSpec - as above for individual lines

if nargin < 1 || isempty(f), f = gcf; end
if nargin < 2 || isempty(axesSpec)
    
    axesSpec = {'FontSize',             14      , ...
                'box',                  'on'    , ...
                'XMinorTick',           'on'    , ...
                'YMinorTick',           'on'    , ...
                'TickLabelInterpreter', 'latex'};
end
if nargin < 3 || isempty(axesSpec)
    
    lineSpec = {'LineWidth',            1       , ...
                'MarkerSize',           14};
end

legendSpec = {'Interpreter',      'latex'};
axes = get(f, 'Children');

for i = 1:numel(axes)
    
    if strcmp(get(axes(i),'tag'), 'legend')
        
        set(axes(i), legendSpec{:});
    else
        set(axes(i), axesSpec{:})
        
        data = get(axes(i), 'Children');
        for j = 1:numel(data)
            
            set(data(j), lineSpec{:})
        end
    end
end

end
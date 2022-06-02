function f = plot_pareto_3to2(cost, axis_labels, invert_cost, highlight)
%PLOT_PARETO takes multi-objective input and plots 3D Pareto fronts in 2D 
%with colourbar for third dimension

axesSpec = {'FontSize',     14      , ...
    'box',                  'on'    , ...
    'XMinorTick',           'on'    , ...
    'YMinorTick',           'on'    , ...
    'TickLabelInterpreter', 'latex'};

if nargin < 3 || isempty(invert_cost)
    
    invert_cost = zeros(1, size(cost, 2));
end

if nargin < 4 || isempty(highlight), highlight = 0; end

normCost = cost;
normCost(:, invert_cost) = -normCost(:, invert_cost);
normCost = normalize(normCost, 'range');

f = figure;
f.Position = [0,353,1050,365];

colour = flip(bone(256));
% Removing too light entries
colormap(colour(25:end,:))

for i = 1:3
    
    subplot(1,3,i)
    hold on
    xlabel(axis_labels{1}, 'Interpreter', 'latex')
    ylabel(axis_labels{2}, 'Interpreter', 'latex')
    set(gca, axesSpec{:});
    grid on
    
    scatter(cost(:,1), cost(:,2), 12, cost(:,3), 'filled')
    
    if any(highlight)
        plot(cost(highlight, 1), cost(highlight, 2), ...
            'ko', 'Markersize', 10, 'LineWidth', 1.5)
    end
    
    % caxis([min(cost(:,3)) max(cost(:,3))])
    cb = colorbar('NorthOutside');
    cb.FontSize = 14;
    cb.TickLabelInterpreter = 'latex';
    cb.Label.Interpreter = 'latex';
    cb.Label.String = axis_labels{3};
    cb.Label.FontSize = 14;
    
    cost        = circshift(cost, 1, 2);
    normCost    = circshift(normCost, 1, 2);
    axis_labels  = circshift(axis_labels, 1, 2);
    axis square
    hold off
end
end
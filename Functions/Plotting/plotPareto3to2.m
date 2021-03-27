function f = plotPareto3to2(cost, axisLabels, invert)
% Plot 3D pareto fronts in 2D with colourbar for third dimension

axesSpec = {'FontSize',     14      , ...
    'box',                  'on'    , ...
    'XMinorTick',           'on'    , ...
    'YMinorTick',           'on'    , ...
    'TickLabelInterpreter', 'latex'};

normCost = cost;
normCost(:,invert) = -normCost(:,invert);
normCost = normalize(normCost, 'range');
f = figure;
f.Position = [0,353,1535,365];

colour = flip(bone(256));
% Removing too light entries
colormap(colour(25:end,:))

for i = 1:3
    
    subplot(1,3,i)
    hold on
    xlabel(axisLabels{1}, 'Interpreter', 'latex')
    ylabel(axisLabels{2}, 'Interpreter', 'latex')
    set(gca, axesSpec{:});
    grid on
    
    scatter(cost(:,1), cost(:,2), 10, cost(:,3), 'filled')
    % caxis([min(cost(:,3)) max(cost(:,3))])
    cb = colorbar;
    cb.FontSize = 14;
    cb.TickLabelInterpreter = 'latex';
    cb.Label.Interpreter = 'latex';
    cb.Label.String = axisLabels{3};
    cb.Label.FontSize = 14;
    
    cost        = circshift(cost, 1, 2);
    normCost    = circshift(normCost, 1, 2);
    axisLabels  = circshift(axisLabels, 1, 2);
    hold off
end
end
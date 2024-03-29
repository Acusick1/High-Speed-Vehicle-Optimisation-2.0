function f = plot_pareto(cost, axisLabels, invert)
%PLOT_PARETO takes multi-objective input and plots Pareto front in
%consistent format

axesSpec = {'FontSize',     14      , ...
    'box',                  'off'    , ...
    'XMinorTick',           'on'    , ...
    'YMinorTick',           'on'    , ...
    'ZMinorTick',           'on'    , ...
    'TickLabelInterpreter', 'latex'};

normCost = cost;
normCost(:,invert) = -normCost(:,invert);
normCost = normalize(normCost, 'range');

f = figure;
f.Position = [0, 42, 720, 520];
hold on
xlabel(axisLabels{1}, 'Interpreter', 'latex')
ylabel(axisLabels{2}, 'Interpreter', 'latex')

if size(cost, 2) == 2
    
    plot(cost(:,1), cost(:,2), 'k.')
else
    scatter3(cost(:,1), cost(:,2), cost(:,3), 12, normCost, 'filled')
    zlabel(axisLabels{3}, 'Interpreter', 'latex')
end

grid on
set(gca, axesSpec{:});
view(45, 25)
%grid on
hold off


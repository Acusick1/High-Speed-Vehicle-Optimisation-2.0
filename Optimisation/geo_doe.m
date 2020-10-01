function [foils, test_var, test_pen] = geo_doe(g_fun, opt_var)
%% DOE
%% TODO: PUT IN CLASS
nOpt = opt_var.nOpt;

% Number of designs to test for geometric violation design of experiment
nDes = 10 * nOpt;
% Number of designs to use for performance design of experiment 
nTest = 1 * nOpt;

% Geometric DOE positions
[~, var] = opt_var.init_lhs(nDes);


% Have to define individually to ensure listeners are added to each object
for i = nDes:-1:1

    foils(i,:) = BezierFoil.define;
end

parts = builder(var, opt_var.nVar, opt_var.name, foils);
% Calculate geometric violations (no performance analysis)
for i = nDes:-1:1
    
    g_fun.value = g_fun.fun(parts{i});
    geo_pen(i,:) = g_fun.penalty;
end

ind = (1:nDes)';
geo_pen = max(0, geo_pen);
mean_geo_pen = mean(geo_pen, 2);

% Sorting for min geometric penalty
col_sort = sortrows([ind, mean_geo_pen], 2);
del = col_sort(nTest+1:end, 1);
col_sort = col_sort(1:nTest, :);

delete(foils(del,:));

% Take corresponding positions to test performance
foils = foils(col_sort(:,1), :);
test_var = var(col_sort(:,1), :);
test_pen = col_sort(:,2);
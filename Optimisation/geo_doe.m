function [foils, test_var, test_pen] = geo_doe(g_fun, var_obj, nDoe)
%% DOE
%% TODO: PUT IN CLASS
nOpt = var_obj.nOpt;

% Number of designs to test for geometric violation design of experiment
nDes = 10000 * nOpt;

% Geometric DOE positions
[~, var] = var_obj.init_lhs(nDes);


% Have to define individually to ensure listeners are added to each object
for i = nDes:-1:1

    foils(i,:) = BezierFoil.define;
end

parts = var_obj.builder(var, foils);
foils = parts{1};
% Calculate geometric violations (no performance analysis)
for i = nDes:-1:1
    
    g_fun.value = g_fun.fun_handle(foils(i));
    geo_pen(i,:) = g_fun.penalties;
end

ind = (1:nDes)';
geo_pen = max(0, geo_pen);
mean_geo_pen = mean(geo_pen, 2);

% Constraint tolerance
con_tol = 0;
% Sorting for min geometric penalty
col_sort = sortrows([ind, mean_geo_pen], 2);
nCandidates = find(col_sort(:,2) <= con_tol, 1, 'last');
candidates = col_sort(1:nCandidates,:);

if nCandidates > nDoe
    
    opt_var_id = var(candidates(:,1), var_obj.opt_var_id);
    crowd = distancecrowding(opt_var_id, candidates(:,2),[],'genotype');
    col_sort2 = sortrows([candidates(:,1), crowd], 2, 'descend');
else
    error('Increase constraint tolerance or number of candidates with goemtetric DOE')
end

winners = col_sort2(1:nDoe, 1);

% For handles
% del = col_sort(nTest+1:end, 1);
% delete(foils(del,:));

% Take corresponding positions to test performance
foils = foils(winners, :);
test_var = var(winners, :);
test_pen = col_sort2(1:nDoe, 2);
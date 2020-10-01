%% Initialise part handles
clear

rng(1)
% rng('shuffle')

save_path = fullfile(pwd, 'Results');
save_dir = 'aerofoilopt';
[~,~,i] = find_file(save_dir, save_path);
save_dir = fullfile(save_path, [save_dir num2str(i)]);
mkdir(save_dir)

nPop = 6;
nfun = 1;
maxIt = 500;

% Have to define individually to ensure listeners are added to each object
for i = nPop:-1:1

    foils(i,:) = BezierFoil.define;
end

parts = {foils};

%% Initialise variables
a = [];
for i = numel(parts):-1:1
    
    [~, b] = parts{i}(1).define;
    
    num(i) = numel([b.val]);
    a = [b a];
end

nVar = sum(num);
var_names = repmat("", 1, nVar);

i = 1;
j = 1;
while true
    
    dim = length(a(j).val);
    var_names((0:dim-1) + i) = a(j).name;
    i = i + dim;
    j = j + 1;
    if i >= nVar, break; end
end

build_fun = @(var) builder(var, num, var_names, foils);

states = Flightstate(deg2rad(5), 5, 10000, []);
aero = Aerodynamics(false);

base = Aerofoil.baseline();
f_v_fun = @(f) [f.area f.checks([9:11 8 5:7])];

base_aero = aero.analyse(states, 1, base);

f_v_base = f_v_fun(base);
f_v = Violation([f_v_base(1:5) 0 nan 0.1], [nan(1, 4) 0.15 nan 0 0.7], f_v_fun, [], "Area");
force_coeffs = [base_aero.force_coeffs];

p_v_fun = @(fc) [fc.Cl];
p_v = Violation(p_v_fun(force_coeffs), nan, p_v_fun, 1, "Cl");

getCost = @(fc) max([fc.Cd]);

cost_fun = @(opt_var)run(opt_var, var, pos, build_fun, aero, states, f_v, p_v, getCost);

var_min = [a.min];
var_max = [a.max];

var = OptVariables(var_min, var_max, var_names);

[doe_foils, doe_var, doe_pen] = geo_doe(f_v, var);
Aerofoil.plotter(doe_foils);
save_figs(save_dir);
[lf_cost, ~, lf_pen] = run(doe_var, [], [], doe_foils, aero, states, f_v, p_v, getCost);

[hf_cost, hf_pen, hf_results] = hfrun(doe_foils, states, p_v, getCost);

% Difference in cost and performance penalties between fidelities
% lf has to be reduced to ensure only performance penalties are compared
data = [hf_cost - lf_cost, hf_pen - lf_pen(:, size(hf_pen, 2)+1:end)];

pso = PSO(var_min, var_max, nPop, cost_fun, [], nfun, maxIt);
save(fullfile(save_dir, 'preopt'))

pso = pso.main();

save(fullfile(save_dir, 'final'))

function [cost, vio, all_vio] = run(opt_var, var, pos, build_fun, aero, states, f_v, p_v, get_cost)
%% Cost function
%% TODO: split into cost a violation functions
%% TODO: don't even bring foils in and rely on handle?
% Or bring in all parts and enter them as varargin in builder

%% Apply variables to objects

[nPop, ~] = size(opt_var);

if isempty(var)

    var = opt_var;
else
    var = repmat(var, nPop, 1);
    var(:, pos) = opt_var;
end

% Nonsense to allow parts to be passed directly
try 
    parts = build_fun(var);
catch
    parts = mat2cell(build_fun, ones(nPop, 1));
end

%% Analysis

feasible = true(nPop, 1);
nStates = numel(states);
Aref = ones(nPop, 1);

% for i = 1:nPop
parfor i = 1:nPop
    
    aeros{i} = aero.analyse(states, Aref(i), parts{i,:});
end

for i = nPop:-1:1
    
    for j = nStates:-1:1
        
        fc = [aeros{i}(j,:).force_coeffs];
        f = [aeros{i}(j,:).forces];
        
        if isempty(fc)
            
            feasible(i,:) = false; 
        else
            % Flipping indices for correct vectorisation later
            fc_sum(j,i) = fc.sum;
            f_sum(j,i) = f.sum;
        end
    end
    
    if feasible(i)
        %% TODO: needs changing from cell outside of foil only opt
        f_v.value = f_v.fun(parts{i,:});
        p_v.value = p_v.fun(fc_sum(:,i));
        all_vio(i,:) = [f_v.penalties, p_v.penalties];
        vio(i,:) = Violation.get_penalty(all_vio(i,:));
        cost(i,:) = get_cost(fc_sum(:,i));
    end
end

cost(~feasible, :) = inf;
vio(~feasible, :) = inf;

end

function [cost, vio, results] = hfrun(foils, states, p_v, get_cost)

aero = SU2();
gmsh = Gmsh();

vals = SU2.convert_states(states);

for i = numel(foils):-1:1
    
    gmsh.generate(foils(i).coords);
    copyfile(gmsh.out_path, aero.path);
    results(i,:) = aero.run(vals);
    cost(i,:) = get_cost(results(i,:));
    
    p_v.value = p_v.fun(results(i,:));
    vio(i,:) = p_v.penalties;
end

end
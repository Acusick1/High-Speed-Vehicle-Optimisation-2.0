classdef PSO < GlobalOptimisation
    
    properties
        
        w = 0.8
        wRange = [0.3, 1.5]
        c1 = 1.49
        c2 = 1.49
        c3 = 0.149
        var_size
        par_vel
        best = struct('variables', [], 'cost', [], 'penalty', [], 'penalties', [])
        gBest
        hist
        gBest_method = "combi"
        max_vel
    end
    properties (Dependent)
        
        fitness
        design_view
    end
    
    methods
        function obj = PSO(lb, ub, nPop, cost_fun, vio_fun, nFun, maxIt)
            
            obj.lb = lb(:)';
            obj.ub = ub(:)';
            obj.nVar = numel(lb);
            obj.mut_prob = min(1/obj.nVar, 0.05);
            obj.nPop = nPop;
            obj.cost_fun = cost_fun;
            
            if nargin >= 5 && ~isempty(vio_fun)
                
                if isobject(vio_fun)
                    
                    obj.vio_fun = repelem(vio_fun, obj.nPop, 1);
                else
                    obj.vio_fun = vio_fun;
                end
            end
            if nargin >= 6 && ~isempty(nFun), obj.nFun = nFun; end
            if nargin >= 7 && ~isempty(maxIt), obj.maxIt = maxIt; end
            
            obj.var_size = [obj.nPop, obj.nVar];
            obj.par_vel = zeros(obj.var_size);
            obj.penalty = zeros(obj.nPop, 1);
        end
        function obj = update_cost(obj)
            %% TODO: Hacky
            % Allows combined cost/vio functions to be used (default) along
            % with cost only functions
            try
                [obj.cost, obj.penalty, obj.penalties] = obj.cost_fun(obj.variables);
            catch
                obj.cost = obj.cost_fun(obj.variables);
            end
            
            if isempty(obj.vio_fun)
                
                obj.penalties = zeros(obj.nPop, 1);
            else
                obj.penalties = max(0, obj.vio_fun(obj.variables));
                obj.penalty = sum(obj.penalties, 2);
            end
        end
        function obj = init_best(obj)
            
            fn = fieldnames(obj.best);
            for i = 1:length(fn)
                
                obj.best.(fn{i}) = obj.(fn{i});
            end
        end
        function a = get.fitness(obj)
            
            cost = obj.cost;
            pen = obj.penalty;
            
            if isempty(pen)
                
                a = cost;
            else
                a = cost + max(0, pen);
            end
        end
        function obj = update_position(obj, hood_pos)
            
            vel = obj.par_vel;
            best_pos = obj.best.variables;
            gBest_pos = obj.gBest.variables;
            
            vel = obj.w * vel + ...
                obj.c1 * rand(obj.var_size).*(best_pos - obj.variables) + ...
                obj.c2 * rand(obj.var_size).*(gBest_pos - obj.variables);
            
            if nargin >= 2 && ~isempty(hood_pos)
                
                vel = vel + ...
                    obj.c3 * rand(obj.var_size).*(hood_pos - obj.variables);
            end
            
            if ~isempty(obj.max_vel)
                
                con1 = vel > obj.max_vel;
                con2 = vel < -obj.max_vel;
                
                vel(con1) = obj.max_vel(con1);
                vel(con2) = -obj.max_vel(con2);
            end
            
            obj.variables = obj.variables + vel;
            obj.par_vel = vel;
            obj = obj.bound();
        end
        function obj = update_best(obj, con_tol)
            
            if nargin < 2 || isempty(con_tol), con_tol = inf; end
            
            id = Population.dominance(obj.cost, obj.best.cost, obj.penalty, obj.best.penalty, con_tol);
            con = id == 1;
            
            fn = fieldnames(obj.best);
            for i = 1:length(fn)
                
                if ~isempty(obj.(fn{i}))
                    
                    obj.best.(fn{i})(con,:) = obj.(fn{i})(con,:);
                end
            end
        end
        function obj = update_gBest(obj, e)
            
            if nargin < 2 || isempty(e), e = 0; end            
            
            best_id = obj.tournament(obj.best.cost, obj.best.penalty, e);
            
            fn = fieldnames(obj.best);
            for i = 1:length(fn)
                
                a.(fn{i}) = obj.best.(fn{i})(best_id,:);
            end
            
            if isempty(obj.gBest) || a.cost < obj.gBest.cost
                
                obj.gBest = a;
            end
        end
        function obj = update_intertia(obj, wc)
             
            % Alter inertia based on above counter
            if wc < 2
                
                obj.w = min(obj.w * 2, max(obj.wRange));
                
            elseif wc > 5
                
                obj.w = max(obj.w * 0.5, min(obj.wRange));
            end
        end
        function obj = bound(obj)
            
            pos = obj.variables;
            
            [lbMat, ubMat] = obj.bMats();
            
            % Ensure new particle position is within boundaries
            con1 = pos > ubMat;
            con2 = pos < lbMat;
            
            % If not set particle velocity to zero and enforce bounds
            obj.par_vel(con1 | con2) = 0;
            pos(con1) = ubMat(con1);
            pos(con2) = lbMat(con2);
            
            obj.variables = pos;
        end
        function a = get.design_view(obj)
            
            [a.best, a.gBest] = ...
                design.from_struct(obj.best, obj.gBest);
        end
        function obj = main_simp(obj)
            
            [obj.variables, obj.best.variables] = ...
                deal(obj.init_lhs(obj.nPop));
            
            obj = obj.update_cost();
            
            obj.best.cost = obj.cost;
            obj.best.penalty = obj.penalty;

            obj = obj.update_gBest();
            obj.hist = obj.gBest;
            
            wc = 0;
            
            for i = 2:obj.maxIt + 1
                
                obj = obj.update_position();
                obj = obj.update_cost();
                obj = obj.update_best();
                obj = obj.update_gBest();
                
                obj.hist(i,:) = obj.gBest;
                
                if any(i-1 == obj.save_it), obj.save_opt(i); end
                
                if isequaln(obj.hist(i), obj.hist(i-1))
                    
                    obj.stall = obj.stall + 1;
                    wc = wc + 1;
                else
                    obj.stall = 0;
                    wc = max(0, wc - 1);
                end
                
                obj = obj.update_intertia(wc);
                
                if obj.verbose
                    
                    fprintf('Iteration %i: f(x): %4.2f p(x): %4.2f\n', i-1, obj.gBest.cost, obj.gBest.penalty);
                end
            end
        end
        function obj = main(obj)
            
            [obj.variables, obj.best.variables] = ...
                deal(obj.init_lhs(obj.nPop));
            
            sets = obj.subset(3);
            
            hood = obj.single_link_neighbour(obj.nPop, 4);
            
            obj = obj.update_cost();
            obj = obj.init_best();
            obj = obj.update_gBest();
            obj = obj.constraint_tolerance();
            
            obj.hist = obj.gBest;
            
            wc = 0;
            
            for i = 2:obj.maxIt + 1
                
                tic
                % Find best fitness/violation value in neighbourhood
                hood_cost = obj.best.cost(hood);
                hood_pen = obj.best.penalty(hood);
                
                id = obj.tournament(hood_cost, hood_pen, obj.con_tol(i));
                
                for j = obj.nPop:-1:1
                    
                    hood_best_var(j,:) = obj.best.variables(hood(j, id(j)), :);
                end
                
                obj = obj.update_position(hood_best_var);
                
                obj.variables(sets(:,2), :) =...
                    obj.uni_mutation(obj.variables(sets(:,2), :));
                
                obj.variables(sets(:,3), :) =...
                    obj.nonuni_mutation(obj.variables(sets(:,3), :), i-1);
                
                obj = obj.update_cost();
                obj = obj.update_best(obj.con_tol(i));
                obj = obj.update_gBest();
                
                obj.hist(i,:) = obj.gBest;
                
                if any(i-1 == obj.save_it), obj.save_opt(i); end
                
                if isequaln(obj.hist(i), obj.hist(i-1))
                    
                    obj.stall = obj.stall + 1;
                    wc = wc + 1;
                else
                    obj.stall = 0;
                    wc = max(0, wc - 1);
                end
                
                obj = obj.update_intertia(wc);
                time = toc;
                
                if obj.verbose
                    
                    fprintf('Iteration %i: f(x): %4.2f p(x): %4.2f e(it): %4.2f t(it): %4.2f\n', i-1, obj.gBest.cost, obj.gBest.penalty, obj.con_tol(i-1), time);
                end
                
                obj.it = i - 1;
                if obj.stall >= obj.maxStall, break; end
            end
        end
    end
    
    methods (Static)
        function obj = test_unconstrained()
            
            nPop = 120;
            cost_fun = @(x) sum(x(:).^2);
            n = 5;
            lb = zeros(n, 1) - 10e9;
            ub = zeros(n, 1) + 10e9;
            obj = PSO(lb, ub, nPop, cost_fun);
            obj = obj.main();
        end
        function obj = test_constrained()
            %% Rosenbrock function constrained with a cubic and a line
            nPop = 120;
            cost_fun = @(x) (1 - x(1)).^2 + 100*(x(2) - x(1).^2).^2;
            vio_fun = @(x) [(x(1)-1).^3 - x(2) + 1, x(1) + x(2) - 2];
            vio_fun = Violation([-1 NaN], [0 0], vio_fun);
            lb = [-1.5 -0.5];
            ub = [1.5 2.5];
            obj = PSO(lb, ub, nPop, cost_fun, vio_fun);
            obj = obj.main_test();
        end
        function hood = single_link_neighbour(nPop, N)
            %% Create singly-linked neighbourhoods
            hood = zeros(nPop, N + 1);
            hood(:,1) = 1:nPop;
            
            for i = 1:nPop
                
                inc = 1;
                neighbours = zeros(1, N);
                
                for j = 1:N
                    
                    if mod(j, 2)
                        
                        neighbour = i + inc;
                    else
                        neighbour = i - inc;
                    end
                    
                    if neighbour > nPop
                        
                        neighbours(j) = neighbour - nPop;
                        
                    elseif neighbour < 1
                        
                        neighbours(j) = nPop + neighbour;
                    else
                        neighbours(j) = neighbour;
                    end
                    
                    inc = inc + 1;
                end
                
                hood(i, 2:end) = neighbours;
            end
        end
    end
end
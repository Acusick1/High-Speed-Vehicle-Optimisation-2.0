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
        function self = PSO(lb, ub, nPop, cost_fun, vio_fun, nFun, max_it)
            %PSO initialisation
            %   Inputs:
            %   lb, ub - lower and upper boundaries of variables
            %   nPop - number of particles in swarm
            %   cost_fun - handle to cost function
            %   vio_fun - handle to violation funciton
            %   nFun - number of cost functions to optimise
            %   max_it - maximum number of iterations
            
            self.lb = lb(:)';
            self.ub = ub(:)';
            self.nVar = numel(lb);
            self.mut_prob = min(1/self.nVar, 0.05);
            self.nPop = nPop;
            self.cost_fun = cost_fun;
            
            %% TODO: Why is this being repeated for every particle? Should be single function like cost_fun
            if nargin >= 5 && ~isempty(vio_fun)
                
                if isobject(vio_fun)
                    
                    self.vio_fun = repelem(vio_fun, self.nPop, 1);
                else
                    self.vio_fun = vio_fun;
                end
            end
            
            if nargin >= 6 && ~isempty(nFun), self.nFun = nFun; end
            if nargin >= 7 && ~isempty(max_it), self.max_it = max_it; end
            
            self.var_size = [self.nPop, self.nVar];
            self.par_vel = zeros(self.var_size);
            self.penalty = zeros(self.nPop, 1);
        end
        
        function self = update_cost(self, i)
            %% TODO: Uncomment when finished testing
            % try
                [self.cost, self.penalty, self.penalties] = self.cost_fun(self.variables);
            % catch ME
            %     obj.save_opt(i);
            %     rethrow(ME)
            % end
            %% TODO: Hacky
            % Allows combined cost/vio functions to be used (default) along
            % with cost only functions
%                 obj.cost = obj.cost_fun(obj.variables);
%                 
%                 if isempty(obj.vio_fun)
%                 
%                     obj.penalties = zeros(obj.nPop, 1);
%                 else
%                     obj.penalties = max(0, obj.vio_fun(obj.variables));
%                     obj.penalty = sum(obj.penalties, 2);
%                 end
%             end
        end
        function self = init_best(self)
            
            fn = fieldnames(self.best);
            for i = 1:length(fn)
                
                self.best.(fn{i}) = self.(fn{i});
            end
        end
        function self = update_position(self, hood_pos)
            
            vel = self.par_vel;
            best_pos = self.best.variables;
            gBest_pos = self.gBest.variables;
            
            vel = self.w * vel + ...
                self.c1 * rand(self.var_size).*(best_pos - self.variables) + ...
                self.c2 * rand(self.var_size).*(gBest_pos - self.variables);
            
            if nargin >= 2 && ~isempty(hood_pos)
                
                vel = vel + ...
                    self.c3 * rand(self.var_size).*(hood_pos - self.variables);
            end
            
            if ~isempty(self.max_vel)
                
                con1 = vel > self.max_vel;
                con2 = vel < -self.max_vel;
                
                vel(con1) = self.max_vel(con1);
                vel(con2) = -self.max_vel(con2);
            end
            
            self.variables = self.variables + vel;
            self.par_vel = vel;
            self = self.bound();
        end
        function self = update_best(self, con_tol)
            
            if nargin < 2 || isempty(con_tol), con_tol = inf; end
            
            id = Population.dominance(self.cost, self.best.cost, self.penalty, self.best.penalty, con_tol);
            con = id == 1;
            
            fn = fieldnames(self.best);
            for i = 1:length(fn)
                
                if ~isempty(self.(fn{i}))
                    
                    self.best.(fn{i})(con,:) = self.(fn{i})(con,:);
                end
            end
        end
        function self = update_gBest(self, e)
            
            if nargin < 2 || isempty(e), e = 0; end            
            
            best_id = self.tournament(self.best.cost, self.best.penalty, e);
            
            fn = fieldnames(self.best);
            for i = 1:length(fn)
                
                a.(fn{i}) = self.best.(fn{i})(best_id,:);
            end
            
            if isempty(self.gBest) || a.cost < self.gBest.cost
                
                self.gBest = a;
            end
        end
        function self = update_intertia(self, wc)
             
            % Alter inertia based on above counter
            if wc < 2
                
                self.w = min(self.w * 2, max(self.wRange));
                
            elseif wc > 5
                
                self.w = max(self.w * 0.5, min(self.wRange));
            end
        end
        function self = bound(self)
            
            pos = self.variables;
            
            [lbMat, ubMat] = self.bMats();
            
            % Ensure new particle position is within boundaries
            con1 = pos > ubMat;
            con2 = pos < lbMat;
            
            % If not set particle velocity to zero and enforce bounds
            self.par_vel(con1 | con2) = 0;
            pos(con1) = ubMat(con1);
            pos(con2) = lbMat(con2);
            
            self.variables = pos;
        end
        function self = run_simple(self)
            
            [self.variables, self.best.variables] = ...
                deal(self.init_lhs(self.nPop));
            
            self = self.update_cost(0);
            
            self.best.cost = self.cost;
            self.best.penalty = self.penalty;

            self = self.update_gBest();
            self.hist = self.gBest;
            
            wc = 0;
            
            for i = 2:self.max_it + 1
                
                self = self.update_position();
                self = self.update_cost(i-1);
                self = self.update_best();
                self = self.update_gBest();
                
                self.hist(i,:) = self.gBest;
                
                if any(i-1 == self.save_it), self.save_opt(i-1); end
                
                if isequaln(self.hist(i), self.hist(i-1))
                    
                    self.stall = self.stall + 1;
                    wc = wc + 1;
                else
                    self.stall = 0;
                    wc = max(0, wc - 1);
                end
                
                self = self.update_intertia(wc);
                
                if self.verbose
                    
                    fprintf('Iteration %i: f(x): %4.2f p(x): %4.2f\n', i-1, self.gBest.cost, self.gBest.penalty);
                end
            end
        end
        function self = run(self)
            %RUN particle swarm optimisation
            
            % Initialise particle variables and particle best variables
            % using Latin hypercube sampling
            self.variables = self.init_lhs(self.nPop);
            self.best.variables = self.variables;
            
            % Initialise mutation subset indices
            sets = self.subset(3);
            
            % Initialise particle neighbourhood indices 
            hood = self.single_link_neighbour(self.nPop, 4);
            
            % Initial PSO loop
            self = self.update_cost(0);
            self = self.init_best();
            self = self.update_gBest();
            self = self.constraint_tolerance();
            
            self.hist = self.gBest;
            
            wc = 0;
            
            for i = 2:self.max_it + 1
                
                tic
                % Find best fitness/violation value in neighbourhood
                hood_cost = self.best.cost(hood);
                hood_pen = self.best.penalty(hood);
                
                id = self.tournament(hood_cost, hood_pen, self.con_tol(i));
                
                for j = self.nPop:-1:1
                    
                    hood_best_var(j,:) = self.best.variables(hood(j, id(j)), :);
                end
                
                self = self.update_position(hood_best_var);
                
                self.variables(sets(:,2), :) =...
                    self.uni_mutation(self.variables(sets(:,2), :));
                
                self.variables(sets(:,3), :) =...
                    self.nonuni_mutation(self.variables(sets(:,3), :), i-1);
                
                self = self.update_cost(i-1);
                self = self.update_best(self.con_tol(i));
                self = self.update_gBest();
                
                self.hist(i,:) = self.gBest;
                
                if any(i-1 == self.save_it)
                
                    self.save_opt(i-1);
                    % Resetting parpool to avoid OOM errors
                    if ~isempty(gcp('nocreate')) && self.max_it >= 500
                        
                        delete(gcp); 
                    end
                end
                
                if isequaln(self.hist(i), self.hist(i-1))
                    
                    self.stall = self.stall + 1;
                    wc = wc + 1;
                else
                    self.stall = 0;
                    wc = max(0, wc - 1);
                end
                
                self = self.update_intertia(wc);
                time = toc;
                
                if self.verbose
                    
                    fprintf('Iteration %i: f(x): %4.2f p(x): %4.2f e(it): %4.2f t(it): %4.2f\n', i-1, self.gBest.cost, self.gBest.penalty, self.con_tol(i-1), time);
                end
                
                self.it = i - 1;
                if self.stall >= self.maxStall, break; end
            end
        end
        function fit = get.fitness(obj)
            
            cost = obj.cost;
            pen = obj.penalty;
            
            if isempty(pen)
                
                fit = cost;
            else
                fit = cost + max(0, pen);
            end
        end
%         function design_view = get.design_view(obj)
%             
%             [design_view.best, design_view.gBest] = ...
%                 Design.from_struct(obj.best, obj.gBest);
%         end
    end
    
    methods (Static)
        function obj = test_unconstrained()
            
            nPop = 120;
            cost_fun = @(x) sum(x(:).^2);
            n = 5;
            lb = zeros(n, 1) - 10e9;
            ub = zeros(n, 1) + 10e9;
            obj = PSO(lb, ub, nPop, cost_fun);
            obj = obj.run();
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
            obj = obj.run();
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
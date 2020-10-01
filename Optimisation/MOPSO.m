classdef MOPSO < PSO
    
    properties
        
        PF = struct('variables', [], 'cost', [], 'penalty', [], 'nDom', [])
        maxPF = 100
    end
    
    methods
        function obj = main(obj)

            obj.variables = deal(obj.init_lhs(obj.nPop));
            
            sets = obj.subset(3);

            obj = obj.update_cost();
            obj = obj.constraint_tolerance();
            
            e = obj.con_tol(1);
            
            obj = obj.update_PF(e);
            
            best_id(1) = obj.roulette(obj.PF.nDom);
            best_id(2) = obj.roulette(obj.PF.nDom, true);
            
            obj = obj.init_best();
            obj = obj.update_gBest(best_id);
            
            % Print iteration, mean for every fitness function, number of PF particles
            str = repmat('%4.3f ', 1, obj.nFun);
            
            obj.hist = obj.gBest;
            
            fn = fieldnames(obj.gBest);
            con = 1:length(best_id);
            
            for j = 1:length(fn)
                
                obj.hist(1).(fn{j})(con,:) = obj.gBest.(fn{j})(con,:);
            end
            
            save(fullfile(obj.save_dir, 'Init'));
            
            for i = 2:obj.maxIt + 1
                
                tic
                obj = obj.update_position();
                
                obj.variables(sets(:,2), :) =...
                    obj.uni_mutation(obj.variables(sets(:,2), :));
                
                obj.variables(sets(:,3), :) =...
                    obj.nonuni_mutation(obj.variables(sets(:,3), :), i);
                
                obj = obj.update_cost();
                
                e = obj.con_tol(i);
                obj = obj.update_PF(e);
                
                best_id(1) = obj.roulette(obj.PF.nDom);
                best_id(2) = obj.roulette(obj.PF.nDom, true);
                
                obj = obj.update_best(e);
                obj = obj.update_gBest(best_id);
                
                for j = 1:length(fn)
                    
                    obj.hist(i,:).(fn{j})(con,:) = obj.gBest.(fn{j})(con,:);
                end
                
                time = toc;
                
                fprintf(['Iteration %i: Mean PF f(x): ' str ' p(x): %4.2f e: %4.2f nPF: %i Success rate: %i t(it): %4.2f\n'], i-1, mean(obj.PF.cost, 1), mean(obj.PF.penalty), obj.conTol(i-1), size(obj.PF.cost, 1), sum(all(isfinite(obj.cost), 2))/obj.nPop*100, time);
                
                if any(i-1 == obj.saveIt)
                    
                    obj.save_opt(i);
                end
            end
        end
        function obj = update_PF(obj, con_tol)

            if nargin < 2 || isempty(con_tol), con_tol = inf; end
            
            %% TODO: This or below
            [nPF, ~] = size(obj.PF.cost);
            external = false(nPF, 1);
            for i = 1:nPF
                
                external(i) = ~any(all(obj.PF.variables(i,:) == obj.best.variables, 2));
            end
            
            fn = fieldnames(obj.best);
            for i = 1:length(fn)
                %% TODO: This or above
                % Only need to compare current population with pareto front
                % Including best leads to unnecessary duplicates
                %% TODO: Ensure PF are saved as bests?
                combi.(fn{i}) = [obj.(fn{i}); obj.best.(fn{i}); ...
                    obj.PF.(fn{i})(external,:)];
            end
            
            feasible = combi.penalty <= con_tol;
            
            [obj.PF.cost, PF_ID, obj.PF.nDom] = ...
                population.nondominated(combi.cost, feasible);
            
            obj.PF.variables = combi.variables(PF_ID,:);
            obj.PF.penalty = combi.penalty(PF_ID,:);
            
            if size(obj.PF.cost, 1) > obj.maxPF
                
                [~, id] = population.maxpf(obj.PF.cost, obj.maxPF);
                
                fn = fieldnames(obj.PF);
                for i = 1:length(fn)
                    
                    obj.PF.(fn{i}) = obj.PF.(fn{i})(id, :);
                end
            end
        end
        function obj = update_gBest(obj, best_id)
            
            if nargin < 2
                %% TODO: Single gBest
                % best_id = pso.roulette(obj.PF.cost, true);
            end
            
            % Replicate to ensure size consistency with population
            rep = obj.nPop/length(best_id);
            
            fn = fieldnames(obj.PF);
            for i = 1:length(fn)
                
                obj.gBest.(fn{i}) = ...
                    repmat(obj.PF.(fn{i})(best_id,:), rep, 1);
            end
        end
    end
end
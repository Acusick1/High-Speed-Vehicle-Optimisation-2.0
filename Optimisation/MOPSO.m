classdef MOPSO < PSO
    
    properties
        
        PF = struct('variables', [], 'cost', [], 'penalty', [], ...
            'penalties', [], 'nDom', [])
        maxPF = 100
    end
    
    methods
        function obj = main(obj)

            obj.variables = deal(obj.init_lhs(obj.nPop));
            
            sets = obj.subset(3);

            obj = obj.update_cost(0);
            obj = obj.constraint_tolerance();
            
            e = obj.con_tol(1,:);
            
            obj = obj.update_PF(e);
            
            best_id(1) = obj.roulette(obj.PF.nDom);
            best_id(2) = obj.roulette(obj.PF.nDom, true);
            
            obj = obj.init_best();
            obj = obj.update_gBest(best_id);
            
            % Print iteration, mean for every fitness function, number of PF particles
            str = repmat('%4.3f ', 1, obj.nFun);
            
            obj.hist = obj.gBest;
            
            fn = fieldnames(obj.gBest);
            for j = 1:length(fn)
                
                obj.hist(1).(fn{j}) = obj.PF.(fn{j})(best_id,:);
            end
            
            if ~isempty(obj.save_it)
                
                save(fullfile(obj.save_dir, 'Init'));
            end
            
            for i = 2:obj.maxIt + 1
                
                tic
                obj = obj.update_position();
                
                obj.variables(sets(:,2), :) =...
                    obj.uni_mutation(obj.variables(sets(:,2), :));
                
                obj.variables(sets(:,3), :) =...
                    obj.nonuni_mutation(obj.variables(sets(:,3), :), i-1);
                
                obj = obj.update_cost(i-1);
                
                e = obj.con_tol(i,:);
                obj = obj.update_PF(e);
                
                % distance = Population.crowding(obj.PF.cost);
                distance = distancecrowding([],obj.PF.cost);
                best_id(1) = obj.roulette(obj.PF.nDom);
                best_id(2) = obj.roulette(distance, true);
                % best_id(2) = obj.roulette(obj.PF.nDom, true);
                
                obj = obj.update_best(e);
                obj = obj.update_gBest(best_id);
                
                for j = 1:length(fn)
                    
                    obj.hist(i,:).(fn{j}) = obj.PF.(fn{j})(best_id,:);
                end
                
                time = toc;
                
                pf_pen = obj.PF.penalty;
                pf_pen(isnan(pf_pen)) = 0;
                
                if obj.verbose
                    
                    fprintf(['Iteration %i: Mean PF f(x): ' str ' p(x): %4.2f e: %4.2f nPF: %i Success rate: %i t(it): %4.2f\n'], i-1, mean(obj.PF.cost, 1), mean(pf_pen(:)), mean(obj.con_tol(i-1,:)), size(obj.PF.cost, 1), sum(all(isfinite(obj.cost), 2))/obj.nPop*100, time);
                end
                
                if any(i-1 == obj.save_it)
                    
                    obj.save_opt(i-1);
                    % Resetting parpool to avoid OOM errors
                    if ~isempty(gcp('nocreate')) && obj.maxIt >= 500
                        
                        delete(gcp); 
                    end
                end
            end
        end
        function obj = main_loop(obj, it1, nIter, rngState)
            
            if nargin < 2 || isempty(it1), it1 = 2; end
            if nargin >= 3 && ~isempty(nIter)
                
                % If extending simulation, set constraint tolerance to
                % original final value
                obj.maxIt = nIter;
                obj.con_tol(end:nIter+1) = obj.con_tol(end);
            end
            if nargin >= 4 && ~isempty(rngState)
                
                rng(rngState);
            end
            
            fn = fieldnames(obj.gBest);
            sets = obj.subset(3);
            
            str = repmat('%4.3f ', 1, obj.nFun);
            
            for i = it1:obj.maxIt + 1
                
                tic
                obj = obj.update_position();
                
                obj.variables(sets(:,2), :) =...
                    obj.uni_mutation(obj.variables(sets(:,2), :));
                
                obj.variables(sets(:,3), :) =...
                    obj.nonuni_mutation(obj.variables(sets(:,3), :), i-1);
                
                obj = obj.update_cost(i-1);
                
                e = obj.con_tol(i,:);
                obj = obj.update_PF(e);
                
                % distance = Population.crowding(obj.PF.cost);
                % distance = distancecrowding([],obj.PF.cost);
                best_id(1) = obj.roulette(obj.PF.nDom);
                % best_id(2) = obj.roulette(distance, true);
                best_id(2) = obj.roulette(obj.PF.nDom, true);
                
                obj = obj.update_best(e);
                obj = obj.update_gBest(best_id);
                
                for j = 1:length(fn)
                    
                    obj.hist(i,:).(fn{j}) = obj.PF.(fn{j})(best_id,:);
                end
                
                time = toc;
                
                pf_pen = obj.PF.penalty;
                pf_pen(isnan(pf_pen)) = 0;
                
                if obj.verbose
                    
                    fprintf(['Iteration %i: Mean PF f(x): ' str ' p(x): %4.2f e: %4.2f nPF: %i Success rate: %i t(it): %4.2f\n'], i-1, mean(obj.PF.cost, 1), mean(pf_pen(:)), mean(obj.con_tol(i-1,:)), size(obj.PF.cost, 1), sum(all(isfinite(obj.cost), 2))/obj.nPop*100, time);
                end
                
                if any(i-1 == obj.save_it)
                    
                    obj.save_opt(i-1);
                    %% Attempt to reset parpool to avoid OOM errors
                    if ~isempty(gcp('nocreate')) && obj.maxIt >= 500
                        
                        delete(gcp); 
                    end
                end
            end
        end
        function obj = update_PF(obj, con_tol)

            if nargin < 2 || isempty(con_tol), con_tol = inf; end
            
            %% TODO: Ensure PF are saved as bests?
            % Which PF particles are external to best, only include those
            % in combi population to avoid duplicates
            [nPF, ~] = size(obj.PF.cost);
            external = false(nPF, 1);
            for i = 1:nPF
                
                external(i) = ~any(all(obj.PF.variables(i,:) == obj.best.variables, 2));
            end
            
            fn = fieldnames(obj.best);
            for i = 1:length(fn)
            
                combi.(fn{i}) = [obj.(fn{i}); obj.best.(fn{i}); ...
                    obj.PF.(fn{i})(external,:)];
            end
            
            %% TODO: If con_tol == 0, one group?
            pen_is_nan = isnan(combi.penalty);
            groups = unique(pen_is_nan, 'rows');
            
            for i = size(groups, 1):-1:1
                
                con = find(all(pen_is_nan == groups(i,:), 2));
                
                [~, pf_ID_group, pf_nDom_group] = ...
                    Population.nondominated(combi.cost(con,:), combi.penalty(con,:), con_tol);
                
                pf_ID{i} = con(pf_ID_group)';
                pf_nDom{i} = pf_nDom_group';
            end
            
            
            [pf_ID, idx] = sort([pf_ID{:}]);
            pf_nDom = [pf_nDom{:}];
            obj.PF.nDom = pf_nDom(idx)';
            
            fn = fieldnames(combi);
            
            for i = 1:numel(fn)
                
                if ~isempty(combi.(fn{i}))
                    
                    obj.PF.(fn{i}) = combi.(fn{i})(pf_ID,:);
                end
            end
            
            if size(obj.PF.cost, 1) > obj.maxPF
                
                [~, id] = Population.limit_pf(obj.PF.cost, obj.maxPF);
                
                fn = fieldnames(obj.PF);
                for i = 1:length(fn)
                    
                    if ~isempty(obj.PF.(fn{i}))
                        
                        obj.PF.(fn{i}) = obj.PF.(fn{i})(id, :);
                    end
                end
            end
        end
        function obj = update_gBest(obj, best_id)
            
            if nargin < 2
                
                best_id = GlobalOptimisation.roulette(obj.PF.nDom, true);
            end
            
            % Replicate to ensure size consistency with population
            rep = obj.nPop/length(best_id);
            
            fn = fieldnames(obj.PF);
            for i = 1:length(fn)
                
                if ~isempty(obj.PF.(fn{i}))
                    
                    obj.gBest.(fn{i}) = ...
                        repmat(obj.PF.(fn{i})(best_id,:), rep, 1);
                end
            end
        end
    end
    
    methods (Static)
        
        function obj = restart(opt_dir, maxIt)
            %% TODO: Opt_file should be opt folder, otherwise Init overwritten
            %% TODO: Put i into obj, take from there
            % load(fullfile(opt_dir, 'Init'), 'sets', 'fn')
            checkpoints = dir(fullfile(opt_dir, 'it*'));

            for i = numel(checkpoints):-1:1
            
                temp = regexp(checkpoints(i).name, '\d*', 'Match');
                it(i) = str2double(temp{:});
            end
            
            % Plus 1 as PSO loop starts at 2
            % Could be plus 2 to start from next it, but repeating last it
            % incase error occured
            [it, id] = max(it+1);
            load(fullfile(checkpoints(id).folder, checkpoints(id).name), 'obj')
            
            try
                load(fullfile(checkpoints(id).folder, checkpoints(id).name), 'rngState')
                rng(rngState);
            catch
                rngState = [];
            end
            
            if nargin < 2, maxIt = []; end
            
            obj.save_dir = opt_dir;
            obj = obj.main_loop(it, maxIt, rngState);
        end
    end
end
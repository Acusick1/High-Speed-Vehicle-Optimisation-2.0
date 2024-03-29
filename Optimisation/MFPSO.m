classdef MFPSO < PSO
    
    properties
        
        doe
        models
        mf_id
        hf_cost
        hf_penalties
        hf_cost_fun
        hf_completed
        hf_confidence
        data
        it_hf
        hf_min_dist = 10
        hf_nInt = 100
    end
    
    methods
        function obj = MFPSO(hf_cost_fun, varargin)
            
            obj = obj@PSO(varargin{:});
            obj.hf_cost_fun = hf_cost_fun;
        end
        function obj = adjust_cost(obj)
            
            for i = numel(obj.models):-1:1
                
                [pred_data(:,i), pred_ci(:,i)] = ...
                    predict(obj.models{i}, obj.variables);
            end
            
            id = obj.mf_id;
            
            if obj.it_hf > obj.hf_min_dist
                
                iInt = min(obj.it_hf - obj.hf_min_dist, obj.hf_nInt);
                pred = [pred_data(:,1) pred_ci(:,1)];
                pred_data(:,1) = interp1([1 obj.hf_nInt], pred, iInt);
            end
            
            % obj.lf_cost = obj.cost;
            pred_cost = obj.cost + pred_data(:,1);
            obj.cost = pred_cost;
            obj.hf_confidence = pred_ci;
            
            % obj.lf_penalties = obj.penalties(id);
            pred_penalties = obj.penalties(:,id) + pred_data(:,2:end);
            obj.penalties(:,id) = pred_penalties;
            obj.penalty = Violation.get_penalty(obj.penalties, 2);
        end
        function obj = update_gBest(obj, e)
            
            if nargin < 2 || isempty(e), e = 0; end            
            
            best_id = obj.tournament(obj.best.cost, obj.best.penalty, e);
            
            fn = fieldnames(obj.best);
            for i = 1:length(fn)
                
                a.(fn{i}) = obj.best.(fn{i})(best_id,:);
            end
            
            if isempty(obj.gBest) || a.cost < obj.gBest.cost
                
                % Has to be run through hf cost function and rechecked
                % against gBest
                if ~ismember(a.variables, obj.hf_completed, 'rows')
                    %% TODO: Tidy
                    % Also, running low fidelity again is a bit hacky
                    [lf_cost, ~, lf_pen] = obj.cost_fun(a.variables);
                    [a.cost, hf_pen] = obj.hf_cost_fun(a.variables);                    
                    a.penalty = Violation.get_penalty([lf_pen(1:obj.mf_id-1) hf_pen]);
                    
                    if all(isfinite(a.cost))

                        obj.hf_completed(end+1,:) = a.variables;
                        obj.hf_cost(end+1,:) = a.cost;
                        obj.hf_penalties(end+1,:) = hf_pen;
                        obj.data(end+1,:) = [a.cost - lf_cost, hf_pen - lf_pen(obj.mf_id)];
                        obj = obj.update_models();
                    end
                end
                
                % Recheck
                if isempty(obj.gBest) || a.cost < obj.gBest.cost
                    
                    obj.gBest = a;
                end
                
                obj.it_hf = 0;
            else
                obj.it_hf = obj.it_hf + 1;
            end
        end
        function obj = update_models(obj)
        
            for i = size(obj.data, 2):-1:1
                
                obj.models{i} = fitrgp(obj.hf_completed, obj.data(:, i),...
                    'KernelFunction','ardsquaredexponential',...
                    'OptimizeHyperparameters','auto',...
                    'HyperparameterOptimizationOptions',...
                        struct('AcquisitionFunctionName',...
                        'expected-improvement-plus',...
                        'UseParallel',1,'Verbose',0,...
                        'MaxObjectiveEvaluations',100));
            end
            
            close all
        end
        function obj = main(obj)
            
            [obj.variables, obj.best.variables] = ...
                deal(obj.init_lhs(obj.nPop));
            
            sets = obj.subset(3);
            
            hood = obj.single_link_neighbour(obj.nPop, 4);
            
            obj = obj.update_cost(0);
            obj = obj.adjust_cost();
            obj = obj.init_best();
            obj = obj.update_gBest();
            obj = obj.constraint_tolerance();
            
            obj.hist = obj.gBest;
            
            wc = 0;
            
            for i = 2:obj.max_it + 1
                
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
                    obj.nonuni_mutation(obj.variables(sets(:,3), :), i);
                
                obj = obj.update_cost(i-1);
                obj = obj.adjust_cost();
                obj = obj.update_best(obj.con_tol(i));
                obj = obj.update_gBest();
                
                obj.hist(i,:) = obj.gBest;
                
                if any(i-1 == obj.save_it)
                
                    obj.save_opt(i-1);
                    % Resetting parpool to avoid OOM errors
                    if ~isempty(gcp('nocreate')) && obj.max_it >= 500
                        
                        delete(gcp); 
                    end
                end
                
                if isequaln(obj.hist(i), obj.hist(i-1))
                    
                    obj.stall = obj.stall + 1;
                    wc = wc + 1;
                else
                    obj.stall = 0;
                    wc = max(0, wc - 1);
                end
                
                obj = obj.update_intertia(wc);
                time = toc;
                
                fprintf('Iteration %i: f(x): %4.2f p(x): %4.2f e(it): %4.2f t(it): %4.2f\n', i-1, obj.gBest.cost, obj.gBest.penalty, obj.con_tol(i-1), time);
            end
        end
    end
end
classdef DesignOfExperiments
% Design of experiments
    
    properties
        
        fun
        options
        var_min
        var_max
        nVar
        nTest
        positions
        designs
        low_fidelity
        high_fidelity
        models
    end
    
    properties (Constant)
        
        nDes = 10 * DesignOfExperiments.nVar
    end
    
    methods
        function obj = DesignOfExperiments(var_min, var_max, fun, pos)
            
            nVar = length(var_min);
            
            if nargin >= 6
                
                obj.positions = pos;
            else
                obj.positions = var_min' + ...
                    lhsdesign(nDes, nVar, 'criterion', 'maximin') .*...
                    (var_max - var_min)';
            end
            
            obj.nVar = nVar;
            obj.var_min = var_min;
            obj.var_max = var_max;
            obj.fun = fun;
        end
        function obj = analyse(obj)
            
            for i = 1:length(obj.fun)
                
                [cost, pen, all_pen] = feval(obj.fun(i), obj.positions);
            end
        end
        function [obj] = penalty_reduction(obj, keep)
            
            id = (1:obj.nDes)';
            pen = max(0, obj.designs.penalty);
            mean_pen = mean(pen, 2);
            
            % Sorting for min penalty
            col_sort = sortrows([id, mean_pen], 2);
            
            % Take corresponding indices and penalties
            rm = col_sort(keep + 1:end, 1);
            
            p = properties(obj.designs);
            for i = 1:length(p)
                
                obj.designs.(p{i})(rm, :) = [];
            end
                        
            obj.positions(rm, :) = [];
        end
        function obj = create_models(obj, which)
            
            [~, nFid] = size(obj.designs);
            
            for i = nFid:-1:1
                
                cost{i} = obj.designs(:, i).cost;
                pen{i} = max(0, obj.designs(:, i).penalty);
                
                if nargin >= 2
                    
                    pen{i} = pen{i}(:, which);
                end
            end
            
            if nFid == 1
                
                data = [cost, pen];
            
            elseif nFid == 2
                
                data = [cost{2} - cost{1}, pen{2} - pen{1}];
            else
                error("No setup for %i fidelities", nFid)
            end
            
            con = any(isinf(data) | isnan(data), 2);
            obj.positions(con,:) = [];
            obj.designs(con,:) =  [];
            data(con,:) = [];
            
            for i = size(data, 2):-1:1
                
                model{i} = fitrgp(obj.positions, data(:, i),...
                    'KernelFunction','ardsquaredexponential',...
                    'OptimizeHyperparameters','auto',...
                    'HyperparameterOptimizationOptions',...
                        struct('AcquisitionFunctionName',...
                        'expected-improvement-plus',...
                        'UseParallel',1,'Verbose',1));
                
                pred_cost = resubPredict(model{i});
                figure();
                hold on
                plot(data(:, i),'r.');
                plot(pred_cost,'b');
                xlabel('x');
                ylabel('y');
                legend({'data','predictions'}, 'Location','Best');
            end
            
            obj.models = model;
        end
    end
end
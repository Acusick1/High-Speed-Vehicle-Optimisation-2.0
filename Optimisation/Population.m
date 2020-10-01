classdef Population < Design
    
    properties
        
        members
        Size
    end
    
    properties (Constant)
        
        maxPF = 100
    end
    
    methods
        function obj = Population(varargin)
            
            % Ensure each member is a design class instance
            if nargin > 1 || ~isobject(varargin{:})
                
                designs = Design(varargin{:});
            else
                designs = varargin{:};
            end
            
            fn = fieldnames(designs);
            for k = 1:numel(fn)
                
                dim = length(designs(1).(fn{k}));
                obj.(fn{k}) = reshape([designs.(fn{k})], [], dim);
            end
            
            obj.members = designs;
            obj.Size = length(designs);
        end
        function a = update(obj)
            
            if size(obj.cost, 2) == 1
                
                [~, a] = min(obj.cost);
            else
                a = obj.nondominated();
            end
        end
        function sets = subset(obj, nSub, interleave)
            
            if mod(obj.Size, nSub)
                
                error("Subsets not evenly divided")
            end
            
            a = 1:obj.Size;
            
            if nargin >= 3 && interleave
                
                sets = reshape(a', nSub, [])';
            else
                sets = reshape(a, [], nSub);
            end
        end
        function [fitness, violation, pos] = tournament(obj, e)
            
            if nargin < 2, e = inf; end
            
            fitness = obj.cost;
            violation = obj.penalty;
            
            pos = 1:numel(fitness);
            
            if isempty(violation)
                
                contenders = true(obj.members, 1);
            else
                contenders = violation <= e;
            end
            
            if any(contenders)
                % Select lowest fitness value
                fitness = fitness(contenders);
                violation = violation(contenders);
                pos = pos(contenders);
                
                [fitness, id] = min(fitness);
                violation = violation(id);
                pos = pos(id);
            else
                % Select lowest violation
                [violation, id] = min(violation);
                fitness = fitness(id);
                pos = pos(id);
            end
        end
    end
    
    methods (Static)
        
        function [nondom_fit, nondom_ID, pars_dominated] = nondominated(fitness, feasible)
            
            [nPop, ~] = size(fitness);
            
            if nargin < 2 || ~any(feasible), feasible = true(nPop ,1); end
            
            for i = nPop:-1:1
                
                contender = fitness(i,:);
                
                if any(contender == inf) || ~feasible(i), continue; end
                
                % If particle dominated by any other
                con1 = all(fitness <= contender, 2);
                con2 = any(fitness(feasible & con1, :) < contender);
                
                % If no to above, particle is non-dominated
                if ~any(con2)
                    % How many particles does non-dominated particle dominate
                    con1 = all(contender <= fitness, 2);
                    con2 = any(contender < fitness(con1, :), 2);
                    nondom_fit(i,:) = contender;
                    nondom_ID(i,:) = i;
                    pars_dominated(i,:) = sum(con2);
                end
            end
            
            % Remove empty slots
            con = nondom_ID == 0;
            nondom_fit(con, :) = [];
            nondom_ID(con) = [];
            pars_dominated(con) = [];
        end
        function [fitness, id] = limit_pf(fitness, maxPF)
            
            num = (1:size(fitness, 1))';
            crowd = Population.crowding(fitness);
            
            col_sort = sortrows([num, crowd], 2, 'descend');
            fitness = fitness(col_sort(1:maxPF, 1), :);
            id = col_sort(1:maxPF, 1);
        end
        function distance = crowding(fitness, maxf, minf)
            %% Crowding distance calculator
            % Computes how close Pareto Front particles are to other particles on the
            % PF. Used to delete PF particles if number of particles within PF >
            % assigned PFmax
            
            if nargin < 3
                % Max and min non-dominated particle fitness functions
                maxf = max(fitness,[],1);
                minf = min(fitness,[],1);
            end
            
            [nPop, nFun] = size(fitness);
            num = (1:nPop)';
            distance = zeros(nPop, 1);
            
            for i = 1:nFun
                
                fsort = sortrows([num, fitness(:,i)],2);
                distance(fsort(1,1)) = inf;
                distance(fsort(end,1)) = inf;
                
                for j = 2:nPop-1
                    
                    distance(fsort(j), 1) = distance(fsort(j), 1) + ...
                        (fsort(j+1, 2) - fsort(j-1, 2))/(maxf(i) - minf(i));
                end
            end
        end
        function [id] = dominance(f1, f2, v1, v2, e)
            
            [n, ~] = size(f1);
            
            for i = n:-1:1
                
                if  nargin == 2 || all([v1(i) v2(i)] <= e) || v1(i) == v2(i)
                    
                    id(i) = which(f1(i,:), f2(i,:));
                else
                    id(i) = which(v1(i,:), v2(i,:));
                end
            end
            
            function [id] = which(f1, f2)
                
                if all(f1 <= f2, 2) && any(f1 < f2, 2)
                    
                    id = 1;
                    
                elseif all(f2 <= f1, 2) && any(f2 < f1, 2)
                    
                    id = 2;
                else
                    id = randi(2);
                end
            end
        end
    end
end
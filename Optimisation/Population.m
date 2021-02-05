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
        function [fit, pen, pos] = tournament(obj, e)
            
            if nargin < 2, e = inf; end
            
            fit = obj.cost;
            pen = obj.penalty;
            
            pos = 1:numel(fit);
            
            if isempty(pen)
                
                contenders = true(obj.members, 1);
            else
                contenders = all(pen <= e, 2);
            end
            
            if any(contenders)
                % Select lowest fitness value
                fit = fit(contenders);
                pen = pen(contenders);
                pos = pos(contenders);
                
                [fit, id] = min(fit);
                pen = pen(id);
                pos = pos(id);
            else
                % Select lowest violation
                [pen, id] = min(sum(pen, 2));
                fit = fit(id);
                pos = pos(id);
            end
        end
    end
    
    methods (Static)
        
        function [nondom_fit, nondom_ID, pars_dominated] = nondominated(fit, pen, e)
            
            [nPop, ~] = size(fit);
            
            if nargin < 2
                
                feasible = true(nPop ,1); 
            else
                feasible = all(pen <= e | isnan(pen), 2);
            end
            
            if any(feasible)
                
                for i = nPop:-1:1
                        
                    if ~feasible(i), continue; end
                    
                    contender = fit(i,:);
                                        
                    % If particle dominated by any other feasible particle
                    con1 = all(fit <= contender, 2);
                    con2 = any(fit(feasible & con1, :) < contender, 2);
                    
                    % If no to above, particle is non-dominated
                    if ~any(con2)
                        % How many particles does non-dominated particle dominate
                        con1 = all(contender <= fit, 2);
                        con2 = any(contender < fit(con1, :), 2);
                        nondom_fit(i,:) = fit(i,:);
                        nondom_ID(i,:) = i;
                        pars_dominated(i,:) = sum(con2);
                    end
                end
                
                 % Remove empty slots
                con = nondom_ID == 0;
                nondom_fit(con, :) = [];
                nondom_ID(con) = [];
                pars_dominated(con) = [];
            else
                pen(isnan(pen)) = 0;
                [~, nondom_ID] = min(sum(pen, 2));
                nondom_fit = fit(nondom_ID, :);
                pars_dominated = nPop;
            end
        end
        function [fitness, id] = limit_pf(fitness, maxPF)
            
            num = (1:size(fitness, 1))';
            % crowd = Population.crowding(fitness);
            crowd = distancecrowding([],fitness);
            
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
            
            comp = ~(isnan(v1) | isnan(v2));
            vBest = min([sum(v1, 2); sum(v2, 2)]);
            
            [n, ~] = size(f1);
            
            for i = n:-1:1
                
                v1i = v1(i, comp(i,:));
                v2i = v2(i, comp(i,:));
                ei = e(comp(i,:));
                % if both violations are less the e-tolerance, dominance
                % based on cost function, unless one of the violations is
                % the best observed. If both or none equal the best, back
                % to cost function based.
                isBest = any([sum(v1i); sum(v2i)] == vBest, 2);
                fbased = all(v1i <= ei) && all(v2i <= ei) && sum(isBest) ~= 1;
                special = isequal(v1i, v2i);
                
                if nargin == 2 || fbased || special
                    
                    id(i) = which(f1(i,:), f2(i,:));
                else
                    [~,id(i)] = min([sum(v1i), sum(v2i)]);
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
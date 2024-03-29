classdef GlobalOptimisation < Optimisation
    
    properties
        
        nPop
        stall = 0
        maxStall = inf
        mut_prob
        con_tol
    end
    
    methods
        function sets = subset(obj, nSub, interleave)
            
            if mod(obj.nPop, nSub)
                
                error("Subsets not evenly divided")
            end
            
            a = 1:obj.nPop;
            
            if nargin >= 3 && interleave
                
                sets = reshape(a', nSub, [])';
            else
                sets = reshape(a, [], nSub);
            end
        end
        function a = rdm(obj, rows, cols)
            
            if nargin < 2, rows = obj.nPop; end
            if nargin < 3, cols = obj.nVar; end
            
            a = rand(rows, cols);
        end
        function pos = uni_mutation(obj, pos)
            %% Uniform mutation
            % Creates random number matrix same size as variable matrix. Any random
            % number less than mutation probability, the corresponding variable is
            % assigned randomly between min and max bounds
            
            rows = size(pos, 1);
            
            con = obj.rdm(rows) <= obj.mut_prob;
            
            [lbMat, ubMat] = obj.bMats(rows);
            
            pos(con) = unifrnd(lbMat(con), ubMat(con));
        end
        function pos = nonuni_mutation(obj, pos, it)
            %% Non-uniform mutation
            % Timestep based mutation. Large mutations possible at beginning of
            % simulation, reduces to smaller mutations as PSO iterations increase
            
            rows = size(pos, 1);
            
            [lbMat, ubMat] = obj.bMats(rows);
            
            con = obj.rdm(rows) <= obj.mut_prob;
            
            a = obj.rdm(rows);
            r = obj.rdm(rows);
            b = 1;
            
            bool = a < 0.5;
            
            id1 = con & bool;
            id2 = con & ~bool;
            
            gamma = 1 - r.^((1 - it/obj.max_it).^b);
            
            pos(id1) = pos(id1) + ...
                (ubMat(id1) - pos(id1)) .* gamma(id1);
            
            pos(id2) = pos(id2) - ...
                (pos(id2) - lbMat(id2)) .* gamma(id2);
        end
        function [lbMat, ubMat] = bMats(obj, rows)
            
            if nargin < 2 || isempty(rows), rows = obj.nPop; end
            
            lbMat = repmat(obj.lb, rows, 1);
            ubMat = repmat(obj.ub, rows, 1);
        end
        function obj = set.con_tol(obj, val)
            
            if numel(val) == 1
                
                val = repmat(val, obj.max_it + 1, 1);
            end
            
            obj.con_tol = val;
        end
        function obj = constraint_tolerance(obj, tf, ef)
            % Constrained optimization by ε constrained particle swarm optimizer with ε-level control
            
            %% TODO: Why was this here?
            % if ~isempty(obj.con_tol), return; end
            
            % Switching condition limited to 500 iterations > set to ef
            % after
            switchIt = min(ceil(obj.max_it/2), 500);
            
            if nargin < 2 || isempty(tf), tf = switchIt; end
            if nargin < 3 || isempty(ef), ef = 0; end
            
            t = obj.max_it + 1;
            cons = max(obj.penalty(isfinite(obj.penalty)));
            
            e = zeros(t, 1) + inf;
            
            if isempty(cons)
                obj.con_tol = e;
                return
            end
            
            % Needs non-zero value to tend towards
            if ef
                
                tEnd = ef; 
            else
                % Adjusting tf for ef = 0 case. Tend towards non-zero
                % value, then linearly reduce to zero
                tEnd = 1e-3;
                tf = 0.75*tf;
            end
            
            e(1) = 0.5*(mean(cons) + min(cons));
            
            for i = 1:switchIt
                    
                if i > tf
                    % Linear
                    ei = interp1([tf switchIt], [tEnd 0], i);
                    %ei = ef + i * (0 - ef)/tf;
                    
                else
                    % Both from: Constrained Optimization by ε Constrained Differential Evolution with Dynamic ε-Level Control
                    B = log(e(1,:)/tEnd)/tf;
                    ei = e(1,:) .* exp(-B * i);
                    
                    % B = atanh(1 - ef/e0)/tf;
                    % ei = e0 * (1 - tanh(B * t));
                end
                
                e(i + 1, :) = ei;
            end
            
            e(i + 1:end, :) = ef;
            obj.con_tol = e;
        end
    end
    
    methods (Static)
        function id = tournament(fitness, violation, e)
            %% Tournament for selecting best
            % Currently, matrices will be compared row by row to give
            % numel(id) = nRows
            % If vectors, assumed single best is to be found
            
            if nargin < 3, e = 0; end
            
            if isvector(fitness), fitness = fitness(:)'; end
            if isvector(violation), violation = violation(:)'; end
            
            [rows, cols] = size(fitness);
            cols = 1:cols;
            
            contenders = violation <= e;
            
            for i = rows:-1:1
                
                ci = contenders(i,:);
                
                if any(ci)
                    % Select lowest fitness value subject to violation <= e
                    minimise = fitness(i, ci);
                    pos = cols(ci);
                else
                    % Select lowest violation value
                    minimise = violation(i, :);
                    pos = cols;
                end
                
                [~, idi] = min(minimise);
                
                id(i,:) = pos(idi);
            end
        end
        
        function ind = roulette(val, rev)
            %% Roulette wheel global best selection
            % Picks global best with respect to val, higher the value,
            % higher value of being chosen (or lower if rev is true)
            
            dim = numel(val);
            a = (1:dim)';
            
            if nargin == 2 && rev
                
                sorted = sortrows([a, 1./val(:)], 2);
            else
                sorted = sortrows([a, val(:)], 2);
            end
            
            partial = zeros(dim + 1, 1);
            
            for i = 2:dim + 1
                
                partial(i) = sorted(i - 1, 2) + partial(i - 1);
            end
            
            rtotal = rand*partial(end);
            
            if rtotal > 0
                
                ind = sorted(rtotal >= partial(1:end-1) & rtotal <= partial(2:end), 1);
            else
                ind = sorted(rtotal <= partial(1:end-1) & rtotal >= partial(2:end), 1);
            end
            
            ind = ind(1);
        end
        function occurences = test_roulette(val, n, rev)
            
            if nargin < 1 || isempty(val), val = randi(10, 10, 1); end
            if nargin < 2 || isempty(n), n = 1000; end
            if nargin < 3 || isempty(rev), rev = false; end
            
            for i = n:-1:1
                
                ind(i,:) = GlobalOptimisation.roulette(val, rev);
            end
            
            a = 1:length(val);
            occurences = val(:);
            
            for i = a
                
                occurences(i, 2) = sum(ind == i);
            end
            
            occurences = sortrows(occurences, 1);
        end
    end
end
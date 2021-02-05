classdef Optimisation
    
    properties
        
        lb
        ub
        nVar
        nFun = 1
        maxIt = 200
        cost_fun
        vio_fun
        it
        cost
        penalty
        penalties
        variables
        fEval
        problem_spec
        parallel
        save_dir = pwd
        verbose = true
    end
    
    properties (Dependent)
        
        save_it
    end
    
    methods
        function a = init(obj, method, dim)
            
            if nargin < 2 || isempty(method), method = "random"; end
            if nargin < 3 || isempty(dim), dim = 1; end
            
            switch method
                case "random"
                    
                    a = obj.initRandom(dim);
                    
                case "lhs"
                    
                    a = obj.initLhs(dim);
            end
        end
        function a = init_lhs(obj, dim)
            
            a = obj.lb + ...
                lhsdesign(dim, obj.nVar, 'criterion', 'maximin') .* ...
                (obj.ub - obj.lb);
        end
        function a = init_random(obj, dim)
            
            lbound = repmat(obj.lb, dim, 1);
            ubound = repmat(obj.ub, dim, 1);
            a = unifrnd(lbound, ubound, [dim, obj.nVar]);
        end
        function a = get.save_it(obj)
            
            if isempty(obj.save_dir)
                
                a = [];
            else
                a = 0:ceil(obj.maxIt/10):obj.maxIt;
            end
        end
        function save_opt(obj, it)
            
            save(fullfile(obj.save_dir, sprintf('it%i', it-1)));
        end
    end
    
    methods (Static)
        
        function vals = get_uncertainty(val, n, step)
            
            % val: array of values to build uncertainty around
            % n: number of uncertainties at either side of val
            % step: delta to step up/down around val
            
            if nargin < 2 || isempty(n), n = 1; end
            if nargin < 3 || isempty(step), step = 0.05 * val; end
            
            % If only positive steps input, take negative so that
            % uncertainty spans above & below val
            if ~any(step < 0), step = [-step; step]'; end
            
            nVal = numel(val);
            
            % Ensure steps are repeated if each val not given explicit step
            if size(step, 1) < nVal, step = repmat(step, nVal, 1); end
            
            % Create uncertainties with val first, then below, then above
            for i = nVal:-1:1
                
                arr = 1:n;
                back = flip(arr) .* step(i,1);
                forw = arr .* step(i,2);
                        
                vals(i,:) = [0, back, forw] + val(i);
            end
        end
    end
end
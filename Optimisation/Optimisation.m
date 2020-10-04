classdef Optimisation
    
    properties
        
        lb
        ub
        nVar
        nFun = 1
        maxIt = 200
        cost_fun
        vio_fun
        cost
        penalty
        penalties
        variables
        fEval
        problem_spec
        parallel
        save_dir = pwd
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
            
            a = 0:ceil(obj.maxIt/10):obj.maxIt;
        end
        function save_opt(obj, it)
            
            save(fullfile(obj.save_dir, sprintf('it%i', it-1)));
        end
    end
end
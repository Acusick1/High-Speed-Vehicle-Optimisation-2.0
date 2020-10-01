classdef OptVariables < Combinable
    %% Optimisation variables class
    % Inputs: OptVariables(var_min, var_max, name, condition, transform, optimise)
    % var_min & var_max required
    
    properties
        
        var
        var_min
        var_max
        name
        condition
        transform
        optimise
        nOpt
        nVar
        A
        b
        Aeq
        beq
    end
    properties (Dependent)
        
        opt_var
    end
    
    methods
        function obj = OptVariables(varargin)
            
            if nargin >= 1
                
                var_min = varargin{1};
                var_max = varargin{2};
                
                obj.var = (var_min(:) + var_max(:))'/2;
                obj.nVar = length(obj.var);
                constant = var_min == var_max;
                obj.var_min = reshape(var_min(~constant), 1, []);
                obj.var_max = reshape(var_max(~constant), 1, []);
                
                if nargin >= 3, obj.name = varargin{3}(:)'; end
                
                if nargin >= 4 && isstring(varargin{4})
                    
                    obj.condition = varargin{4}(:)';
                else
                    obj.condition = repmat("", 1, obj.nVar);
                end
                
                if nargin >= 5 && isstring(varargin{5})
                    
                    obj.transform = varargin{5}(:)';
                else
                    obj.transform = repmat("", 1, obj.nVar);
                end
                
                if nargin >= 6
                    
                    obj.optimise = varargin{6}(:)';
                else
                    obj.optimise = true(1, obj.nVar);
                end
                
                obj.optimise(constant) = false;
                obj.nOpt = sum(obj.optimise);
            end
        end
        function obj = lincon(obj)
            
            [obj.A, obj.Aeq] = deal(zeros(obj.nVar));
            [obj.b, obj.beq] = deal(zeros(obj.nVar, 1));
            
            for i = 1:obj.nVar
                
                if ~isempty(obj.condition(i))
                    
                    % Eqaulity or inequality (default is inequality)
                    if contains(obj.condition(i), "=") && ~contains(obj.condition(i), "<")
                        
                        [obj.Aeq(i,:), obj.beq(i)] = obj.str2con(i);
                    else
                        [obj.A(i,:), obj.b(i)] = obj.str2con(i);
                    end
                end
            end
            
            if ~any(obj.A(:))
                
                obj.A = [];
                obj.b = [];
            end
            
            if ~any(obj.Aeq(:))
                
                obj.Aeq = [];
                obj.beq = [];
            end
        end
        function [Arow, brow] = str2con(obj, row)
            % Requries equations to be entered with spaces between major statements,
            % number at end must be constraint value
            % Examples: "3 x1 - 3 x2" will assume "<= 0"
            %           "-3 x1 <= 4"
            
            Arow = zeros(1, obj.nVar);
            array = 1:obj.nVar;
            
            splitstring = split(obj.condition(row), ' ');
            
            % Special cases
            if any(contains(splitstring, "previous"))
                
                Arow(row) = 1;
                Arow(row - 1) = -1;
                
            elseif any(contains(splitstring, "next"))
                
                Arow(row) = 1;
                Arow(row + 1) = -1;
            end
            
            % Variable name based
            for i = 1:length(splitstring)
                
                check = str2double(splitstring(i));
                col = array(splitstring(i) == obj.name);
                
                if ~isnan(check)
                    
                    num = check;
                    
                    if i ~= 1 && splitstring(i - 1) == "-"
                        
                        num = -num;
                    end
                    
                elseif ~isempty(col)
                    
                    try
                        Arow(col) = num;
                    catch
                        Arow(col) = 1;
                    end
                end
            end
            
            % Has to be check here
            if ~isnan(check)
                
                brow = check;
            else
                brow = 0;
            end
        end
        function a = get.opt_var(obj)
            
            arr = 1:obj.nVar;
            a = arr(obj.optimise);
        end
        function a = get.nOpt(obj)
            
            a = sum(obj.optimise);
        end
        function a = get.nVar(obj)
            
            a = numel(obj.var);
        end
        function [a, b] = init(obj, method, dim)
            
            if nargin < 2 || isempty(method), method = "random"; end
            if nargin < 3 || isempty(dim), dim = 1; end
            
            switch method
                case "random"
                    
                    [a, b] = obj.init_random(dim);
                    
                case "lhs"
                    
                    [a, b] = obj.init_lhs(dim);
            end
        end
        function [a, b] = init_lhs(obj, dim)
            
            a = [obj.var_min] + ...
                lhsdesign(dim, sum([obj.nOpt]), 'criterion', 'maximin') .* ...
                ([obj.var_max] - [obj.var_min]);
            
            b = obj.combine_var(a);
        end
        function [a, b] = init_random(obj, dim)
            
            lbound = repmat([obj.var_min], dim, 1);
            ubound = repmat([obj.var_max], dim, 1);
            a = unifrnd(lbound, ubound, [dim, sum([obj.nOpt])]);
            
            b = obj.combine_var(a);
        end
        function b = combine_var(obj, opt_var)
            
            b = repmat(obj.var, size(opt_var, 1), 1);
            b(:,obj.opt_var) = opt_var;
        end
    end
    
    methods (Static)
        function a = test_gen()
            
            var_min = [3 -5 0 0];
            var_max = [8 -3 1 1];
            n = ["one", "two", "three","four"];
            c = ["", "", "", "< previous"];
            t = [];
            
            a = OptVariables(var_min, var_max, n, c, t);
            a = a.lincon();
        end
    end
end
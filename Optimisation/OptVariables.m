classdef OptVariables < Combinable
    %% Optimisation variables class
    % Inputs: OptVariables(var_min, var_max, name, condition, transform, optimise)
    % var_min & var_max required
    
    properties
        
        var
        var_min
        var_max
        var_names
        condition
        transform
        optimise
        nOpt
        nVar
        partition
        A
        b
        Aeq
        beq
    end
    properties (Dependent)
        
        opt_var
        view
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
                
                if nargin >= 3, obj.var_names = varargin{3}(:)'; end
                
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
        function parts = builder(obj, var, varargin)
            %% TODO: Varargout to match inputs?
            [nPop, ~] = size(var);
            
            parts = varargin;
            
            n = obj.partition;
            if isempty(n), n = obj.nVar; end
            
            pos_ind = 1:sum(n);
            cumnum = cumsum(n);
            
            next = 1;
            part = 0;
            while true
                
                if part == 0 || next > cumnum(part)
                    
                    part = part + 1;
                end
                
                % Current variable may be vector, so choose all with same
                % variable name 
                current = obj.var_names(next);
                % Variable names may repeat across parts, so only allow
                % variables to be grabbed within current part variable 
                % number bounds
                con = current == obj.var_names;
                
                if part > 1
                    
                    outside = [1:cumnum(part-1), cumnum(part)+1:numel(con)];
                else
                    outside = cumnum(part)+1:numel(con);
                end
                
                con(outside) = false;
                
                vari = mat2cell(var(:,con), ones(nPop, 1), sum(con));
                
                % Hack to ensure single variable set calls can be run
                try
                    [parts{part}.(current)] = vari{:};
                catch
                    parts{part} = parts{part}(1);
                    parts{part}(1).(current) = vari{:};
                end
                
                next = find(~con & pos_ind > next, 1);
                if isempty(next), break; end
            end
            
            %% Change to fully cell based
            % Cannot do earlier due to struct access above
%             for i = numel(parts):-1:1
%                 
%                 new(:,i) = mat2cell(parts{i}, ones(nPop, 1));
%             end
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
                col = array(splitstring(i) == obj.var_names);
                
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
        function a = get.view(obj)
            
            arr = 1:obj.nVar;
            isOpt = ismember(arr, obj.opt_var)';
            arr(isOpt) = 1:obj.nOpt;
            arr(~isOpt) = 0;
            a = table(obj.var_names', arr', 'VariableNames', {'name', 'opt_id'});
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
        function v = sensitivity(obj, id, n)
            %% Sensitivity analysis variable setup
            % Inputs: id (int, int vector)
            %         n  (number of steps)
            vars = numel(id);
            if nargin < 3 || isempty(n), n = 10; end
            
            if numel(n) ~= vars, n = repmat(n, vars, 1); end
            
            dFF = fullfact(n);
            
            minv = obj.var_min;
            maxv = obj.var_max;
            v = repmat((minv + maxv)/2, size(dFF, 1), 1);
            
            for i = 1:vars
                
                j = id(i);
                vsens = minv(j):(maxv(j) - minv(j))/(n(i)-1):maxv(j);
                
                v(:,j) = vsens(dFF(:,i));
            end
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
classdef OptVariables < Combinable
    %% TODO: Split into two classes? One that contains part based optvariables, conditions etc, and one meta-object that takes part objects as an input?
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
        function self = OptVariables(varargin)
            
            if nargin > 0
                
                var_min = varargin{1};
                var_max = varargin{2};
                
                self.var = (var_min(:) + var_max(:))'/2;
                self.nVar = length(self.var);
                constant = var_min == var_max;
                self.var_min = reshape(var_min(~constant), 1, []);
                self.var_max = reshape(var_max(~constant), 1, []);
                
                if nargin >= 3, self.var_names = varargin{3}(:)'; end
                
                if nargin >= 4 && isstring(varargin{4})
                    
                    self.condition = varargin{4}(:)';
                else
                    self.condition = repmat("", 1, self.nVar);
                end
                
                if nargin >= 5 && isstring(varargin{5})
                    
                    self.transform = varargin{5}(:)';
                else
                    self.transform = repmat("", 1, self.nVar);
                end
                
                if nargin >= 6
                    
                    self.optimise = varargin{6}(:)';
                else
                    self.optimise = true(1, self.nVar);
                end
                
                self.optimise(constant) = false;
                self.nOpt = sum(self.optimise);
            end
        end
        function objects = builder(self, vars, varargin)
            %% TODO: Move this to meta-object discuss at top
            %BUILDER applies input variables to objects
            [nObj, ~] = size(vars);
            objects = varargin;
            
            n = cumsum([0 self.partition]);
            if isempty(n), n = self.nVar; end
            
            for i = 1:numel(objects)
                
                object_id = n(i) + 1 : n(i+1);
                object_vars = vars(:, object_id);
                object_var_names = self.var_names(:, object_id);
                
                objects{i} = variables2object(...
                    objects{i}, ...
                    object_vars, ...
                    object_var_names);
            end
            
            function part = variables2object(part, vars, var_names)
            
                unqiue_props = unique(var_names, 'stable');

                for prop_id = 1:numel(unqiue_props)
                    
                    prop = unqiue_props(prop_id);
                    id = prop == var_names;
                    
                    if nObj == 1
                        part.(prop) = vars(:, id);
                    else
                        prop_cell = mat2cell(...
                            vars(:, id), ...
                            ones(nObj, 1), ...
                            sum(id));
                        
                        [part.(prop)] = prop_cell{:};
                    end
                end
            end
        end
        function self = lincon(self)
            
            [self.A, self.Aeq] = deal(zeros(self.nVar));
            [self.b, self.beq] = deal(zeros(self.nVar, 1));
            
            for i = 1:self.nVar
                
                if ~isempty(self.condition(i))
                    
                    % Eqaulity or inequality (default is inequality)
                    if contains(self.condition(i), "=") && ~contains(self.condition(i), "<")
                        
                        [self.Aeq(i,:), self.beq(i)] = self.str2con(i);
                    else
                        [self.A(i,:), self.b(i)] = self.str2con(i);
                    end
                end
            end
            
            if ~any(self.A(:))
                
                self.A = [];
                self.b = [];
            end
            
            if ~any(self.Aeq(:))
                
                self.Aeq = [];
                self.beq = [];
            end
        end
        function [Arow, brow] = str2con(self, row)
            % Requries equations to be entered with spaces between major statements,
            % number at end must be constraint value
            % Examples: "3 x1 - 3 x2" will assume "<= 0"
            %           "-3 x1 <= 4"
            
            Arow = zeros(1, self.nVar);
            array = 1:self.nVar;
            
            splitstring = split(self.condition(row), ' ');
            
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
                col = array(splitstring(i) == self.var_names);
                
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
        function a = get.opt_var(self)
            
            arr = 1:self.nVar;
            a = arr(self.optimise);
        end
        function a = get.nOpt(self)
            
            a = sum(self.optimise);
        end
        function a = get.nVar(self)
            
            a = numel(self.var);
        end
        function a = get.view(self)
            
            arr = 1:self.nVar;
            isOpt = ismember(arr, self.opt_var)';
            arr(isOpt) = 1:self.nOpt;
            arr(~isOpt) = 0;
            a = table(self.var_names', arr', 'VariableNames', {'name', 'opt_id'});
        end
        function [a, b] = init(self, method, dim)
            
            if nargin < 2 || isempty(method), method = "random"; end
            if nargin < 3 || isempty(dim), dim = 1; end
            
            switch method
                case "random"
                    
                    [a, b] = self.init_random(dim);
                    
                case "lhs"
                    
                    [a, b] = self.init_lhs(dim);
            end
        end
        function [a, b] = init_lhs(self, dim)
            
            a = [self.var_min] + ...
                lhsdesign(dim, sum([self.nOpt]), 'criterion', 'maximin') .* ...
                ([self.var_max] - [self.var_min]);
            
            b = self.combine_var(a);
        end
        function [a, b] = init_random(self, dim)
            
            lbound = repmat([self.var_min], dim, 1);
            ubound = repmat([self.var_max], dim, 1);
            a = unifrnd(lbound, ubound, [dim, sum([self.nOpt])]);
            
            b = self.combine_var(a);
        end
        function v = sensitivity(self, id, n)
            %% Sensitivity analysis variable setup
            % Inputs: id (int, int vector)
            %         n  (number of steps)
            vars = numel(id);
            if nargin < 3 || isempty(n), n = 10; end
            
            if numel(n) ~= vars, n = repmat(n, vars, 1); end
            
            dFF = fullfact(n);
            
            minv = self.var_min;
            maxv = self.var_max;
            v = repmat((minv + maxv)/2, size(dFF, 1), 1);
            
            for i = 1:vars
                
                j = id(i);
                vsens = minv(j):(maxv(j) - minv(j))/(n(i)-1):maxv(j);
                
                v(:,j) = vsens(dFF(:,i));
            end
        end
        function b = combine_var(self, opt_var)
            
            b = repmat(self.var, size(opt_var, 1), 1);
            b(:,self.opt_var) = opt_var;
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
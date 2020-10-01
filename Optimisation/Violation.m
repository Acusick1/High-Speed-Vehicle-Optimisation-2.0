classdef Violation < handle
    
    properties (SetObservable)
       
        value
    end
    
    properties
        
        val_min
        val_max
        scale
        description
        fun
        penalty
        penalties
    end

    properties (Dependent)
        
        view
    end
    
    methods
        function obj = Violation(val_min, val_max, fun, scale, description)
            %% Constructor
            if nargin > 0
                
                obj.val_min = val_min(:)';
                obj.val_max = val_max(:)';
                
                n = numel(val_min);
                
                obj.scale = zeros(1, n);
                obj.description = string.empty(0, n);
                
                if nargin >= 3 && ~isempty(fun)
                    
                    obj.fun = fun;
                else
                    obj.fun = @(x) x;
                end
                
                if nargin < 4 || isempty(scale)
                    
                    obj.scale = nan(size(val_min));
                else
                    obj.scale = scale;
                end
                
                if nargin >= 5, obj.description = description; end
                
                obj = obj.get_scale();
            end
            
            addlistener(obj, 'value', 'PostSet', @obj.evaluate);
        end
        function obj = get_scale(obj)
            % Look for scaling parameter if scaling parameter is empty
            temp = obj.scale;
            
            % Check validity
            rescale = temp == 0 | isnan(temp);
            
            % First option is mean maximum/minimum, then maximum boundary,
            % then minimum boundary
            minCon = ~isnan(obj.val_min);
            maxCon = ~isnan(obj.val_max);
            avg = (obj.val_min + obj.val_max)/2;
            
            % Check to see if mean can be used as scaling factor
            con = rescale &...
                isnumeric(avg) & isfinite(avg) & avg ~= 0;
            temp(con) = avg(con);
            rescale(con) = false;
            
            % Else if max can be used
            con = rescale & maxCon & obj.val_max ~= 0;
            temp(con) = obj.val_max(con);
            rescale(con) = false;
            
            % Else if min can be used
            con = rescale & minCon & obj.val_min ~= 0;
            temp(con) = obj.val_min(con);
            rescale(con) = false;
            
            % Else final option > no rescaling
            temp(rescale) = 1;
            
            obj.scale = temp;
        end
        function obj = evaluate(obj, varargin)
            
%             if nargin > 1 || ~isempty(varargin)
%                 
%                 val = obj.fun(varargin{:});
%             else
                val = obj.value;
%             end
            
            mini = obj.val_min;
            maxi = obj.val_max;
            sc = obj.scale;
            
            % Impose conditions
            con(1,:) = mini - val;
            con(2,:) = val - maxi;
            
            % Flip signs of any with negative scales
            flip = sc < 0;
            con(:, flip) = -con(:, flip);
            
            % Normalise wrt scaling value
            all = con./sc;
            all(isnan(all)) = [];

            obj.penalty = obj.get_penalty(all);
            obj.penalties = all(:)';
            % obj.value = val;
        end
        
        function a = get.view(obj)
            
            vio = obj.penalties;
            
            a(length(obj.val_max), :) = struct();
            
            i = 1;
            for j = 1:length(vio)
                
                if mod(j, 2)
                    
                    if ~isempty(obj.description)
                        
                        a(i).description = obj.description(i);
                    end
                    
                    con = obj.val_min(i);
                    
                    if ~isnan(con)
                        
                        a(i).val_minVio = vio(j);
                        a(i).val_min = con;
                    end
                    
                    a(i).value = obj.value(i);
                else
                    con = obj.val_max(i);
                   
                    if ~isnan(con)
                        
                        a(i).val_max = con;
                        a(i).val_maxVio = vio(j);
                    end
                    
                    i = i + 1;
                end
            end
        end
    end
    
    methods (Static)
        function pen = flatten(pen)
            
            pen = max(0, pen);
        end
        function pen = get_penalty(penalties)
            
            active = penalties(penalties > 0);
            
            % Quad loss
            % a = eta * sum(active(:).^2);
            % SCV
            pen = sum(active(:));
        end
        function obj = test()
            
            val_min = [0 -1];
            val_max = [1 0];
            
            obj = Violation(val_min, val_max, @(x, y) [x y]' + rand(2, 1));
            obj.value = [1, -1];
        end
    end
end
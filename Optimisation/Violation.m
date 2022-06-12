classdef Violation < Combinable% < handle
    
    properties %(SetObservable)
       
        value
    end
    
    properties
        
        val_min
        val_max
        scale
        description
        penalty
        penalties
        fun_handle
    end

    properties (Dependent)
        
        view
    end
    
    methods
        function self = Violation(val_min, val_max, scale, description, fun_handle)
            %% Constructor
            if nargin > 0
                
                self.val_min = val_min(:)';
                self.val_max = val_max(:)';
                
                n = numel(val_min);
                
                self.scale = zeros(1, n);
                self.description = string.empty(0, n);
                
                if nargin < 3 || isempty(scale)
                    
                    self.scale = nan(size(val_min));
                else
                    self.scale = scale;
                end
                
                if nargin >= 4, self.description = description; end
                if nargin >= 5, self.fun_handle = fun_handle; end
                
                self = self.get_scale();
            end
            
            % addlistener(obj, 'value', 'PostSet', @obj.evaluate);
        end
        function self = set.value(self, val)
            
            self.value = val;
            self = self.evaluate();
        end
        function self = get_scale(self)
            % Look for scaling parameter if scaling parameter is empty
            temp = self.scale;
            
            % Check validity
            rescale = temp == 0 | isnan(temp);
            
            % First option is mean maximum/minimum, then maximum boundary,
            % then minimum boundary
            minCon = ~isnan(self.val_min);
            maxCon = ~isnan(self.val_max);
            avg = (self.val_min + self.val_max)/2;
            
            % Check to see if mean can be used as scaling factor
            con = rescale &...
                isnumeric(avg) & isfinite(avg) & avg ~= 0;
            temp(con) = avg(con);
            rescale(con) = false;
            
            % Else if max can be used
            con = rescale & maxCon & self.val_max ~= 0;
            temp(con) = self.val_max(con);
            rescale(con) = false;
            
            % Else if min can be used
            con = rescale & minCon & self.val_min ~= 0;
            temp(con) = self.val_min(con);
            rescale(con) = false;
            
            % Else final option > no rescaling
            temp(rescale) = 1;
            
            self.scale = temp;
        end
        function [self, penalty, penalties] = evaluate(self)
            
            val = self.value;
            mini = self.val_min;
            maxi = self.val_max;
            sc = self.scale;
            
            % Impose conditions
            con(1,:) = mini - val;
            con(2,:) = val - maxi;
            
            % Flip signs of any with negative scales
            flip = sc < 0;
            con(:, flip) = -con(:, flip);
            
            % Normalise wrt scaling value
            all = con./sc;
            all(isnan(all)) = [];

            penalty = self.get_penalty(all);
            penalties = all(:)';
            
            self.penalty = penalty;
            self.penalties = penalties;
        end
        function view = get.view(self)
            
            vio = self.penalties;
            
            view(length(self.val_max), :) = struct();
            
            i = 1;
            for j = 1:length(vio)
                
                if mod(j, 2)
                    
                    if ~isempty(self.description)
                        
                        view(i).description = self.description(i);
                    end
                    
                    con = self.val_min(i);
                    
                    if ~isnan(con)
                        
                        view(i).val_minVio = vio(j);
                        view(i).val_min = con;
                    end
                    
                    view(i).value = self.value(i);
                else
                    con = self.val_max(i);
                   
                    if ~isnan(con)
                        
                        view(i).val_max = con;
                        view(i).val_maxVio = vio(j);
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
        function pen = get_penalty(penalties, dim, scale)
            
            if nargin < 2 || isempty(dim), dim = 1; end
            if nargin < 3 || isempty(scale), scale = false; end
            
            if isvector(penalties)
                
                penalties = penalties(:);
            end
            
            active = max(penalties, 0);
            % Quad loss
            % a = eta * sum(active.^2, dim);
            % SCV
            pen = sum(active, dim);
            
            if scale, pen = pen./(size(penalties, dim)); end
        end
        function obj = test()
            
            val_min = [0 -1];
            val_max = [1 0];
            
            obj = Violation(val_min, val_max, @(x, y) [x y]' + rand(2, 1));
            obj.value = [1, -1];
        end
    end
end
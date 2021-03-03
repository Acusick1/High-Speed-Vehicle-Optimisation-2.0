classdef Combinable
    
    properties
        
    end
    
    methods
        
        function obj = combine(objarray)
            
            fn = fieldnames(objarray, '-full');
            
            obj = objarray(1);
            
            for i = 1:numel(fn)
                
                try
                obj.(fn{i}) = [objarray.(fn{i})];
                catch
                end
            end
        end
        function obj = sum(objarray)
            
            fn = fieldnames(objarray, '-full');
            
            obj = objarray(1);
            
            for i = 1:numel(fn)
               
                obj.(fn{i}) = sum([objarray.(fn{i})]);
            end
        end
        function obj = mean_sigma_var(obj, sigma)
            
            if numel(obj) > 1, obj = obj.combine; end
            if nargin < 2 || isempty(sigma), sigma = 3; end
            
            fn = fieldnames(obj);
            
            for i = 1:numel(fn)
                
                val = obj.(fn{i});
                obj.(fn{i}) = mean(val) + sigma*var(val);
            end
        end
    end
    
    methods (Static)
        
        function obj = hori(varargin)
                        
            objarray = [varargin{:}];
            
            fn = fieldnames(objarray(1), '-full');
            
            obj = objarray(1);
            
            for i = 1:numel(fn)
                
                try
                    obj.(fn{i}) = [objarray.(fn{i})];
                catch
                end
            end
        end
        function obj = vert(varargin)
            
            for i = 1:nargin
                
                varargin{i} = structfun(@transpose, varargin{i}, 'UniformOutput', false);
            end
            
            obj = Combinable.hori(varargin{:});
            
            obj = structfun(@transpose, obj, 'UniformOutput', false);
        end
    end
end
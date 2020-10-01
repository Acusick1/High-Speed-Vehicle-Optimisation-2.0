classdef Combinable
    
    properties
        
    end
    
    methods
        
        function obj = combine(objarray)
            
            fn = fieldnames(objarray);
            
            obj = objarray(1);
            
            for i = 1:numel(fn)
               
                obj.(fn{i}) = [objarray.(fn{i})];
            end
        end
        function obj = sum(objarray)
            
            fn = fieldnames(objarray, '-full');
            
            obj = objarray(1);
            
            for i = 1:numel(fn)
               
                obj.(fn{i}) = sum([objarray.(fn{i})]);
            end
        end
    end
end
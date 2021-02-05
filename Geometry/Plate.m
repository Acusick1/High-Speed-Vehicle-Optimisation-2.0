classdef Plate < Body
    
    properties
        
        length = 1
    end
    
    methods
        function obj = Plate()
            
            obj.conical = false;
        end
        function dogenerate(obj)
            
            yDisc = obj.yDisc - 0.5;
            rows = numel(obj.xDisc);
            cols = numel(yDisc);
            x = obj.xDisc * obj.length;
            y = ones(rows, 1) * yDisc;
            z = zeros(rows, cols);
            
            obj.x = repmat(x, 1, cols);
            obj.y = y;
            obj.z = z;
        end
    end
end
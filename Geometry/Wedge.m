classdef Wedge < Body
    
    properties
        
        length = 1
        angle
    end
    
    methods
        function obj = Wedge(angle)
            
            obj.conical = false;
            
            if nargin >= 1
                
                obj.angle = angle;
            end
        end
        function dogenerate(obj)
            
            yDisc = obj.yDisc - 0.5;
            rows = numel(obj.xDisc);
            cols = numel(yDisc);
            x = obj.xDisc * obj.length;
            y = ones(rows, 1) * [yDisc, flip(yDisc)];
            z = x .* [ones(1, cols) * sin(obj.angle), ones(1, cols) * -sin(obj.angle)];
            
            obj.x = repmat(x, 1, cols*2);
            obj.y = y;
            obj.z = z;
        end
    end
end
classdef Cone < Body
    
    properties
        
        length = 1
        angle
    end
    
    methods
        function obj = Cone(angle)
            
            if nargin >= 1
                
                obj.angle = angle;
            end
        end
        function generate(obj)
            
            yDisc = deg2rad(0:2:180);
            xDisc = (0:0.001:1)';
            % xDisc = obj.xDisc;
            x = xDisc * obj.length;
            y = x*sin(obj.angle) .* sin(yDisc);
            z = x*sin(obj.angle) .* cos(yDisc);
            
            obj.x = repmat(x, 1, numel(yDisc));
            obj.y = y;
            obj.z = z;
        end
    end
end
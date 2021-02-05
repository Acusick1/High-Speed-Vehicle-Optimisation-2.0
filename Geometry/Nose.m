classdef Nose < Geometry
    
    properties (SetObservable)
        radius
        length            % Non-dimensional > multiplied by radius
        offset = [0 0]
        rotation = 0
        xPanels = 10
        yPanels = 10
    end
    
    methods
        function obj = Nose(radius, length, offset)
            
            if nargin > 0
                
                obj.radius = radius;
                obj.length = length;
                
                if nargin >= 3
                    
                    obj.offset = offset;
                end
            end
            
            addlistener(obj,'radius','PostSet',@obj.initialise);
            addlistener(obj,'length','PostSet',@obj.initialise);
            addlistener(obj,'offset','PostSet',@obj.initialise);
            addlistener(obj,'rotation','PostSet',@obj.initialise);
            addlistener(obj,'yPanels','PostSet',@obj.redisc);
        end
        function obj = dogenerate(obj)
            
            dTheta = pi/obj.yPanels;
            phi = 0:dTheta:pi;
            
            len = obj.length;
            rad = obj.radius;
            dim = size(phi);
            
            if all([rad, len])
                
                max_theta = acos(-(len/rad)+1);
                theta = (0:max_theta/obj.xPanels:max_theta)';
                
                x = rad * (1 - cos(theta)) * ones(dim);
                y = rad * sin(theta) .* sin(phi);
                z = rad * cos(phi) .* sin(theta);
                
            else
                [x,y,z] = deal(zeros(dim));
            end
            
            obj.x = x;
            obj.y = y;
            obj.z = z;
            
            if any(obj.offset)
                
                % obj.rotate(obj.rotation);
                obj.translate(obj.offset(1), [], obj.offset(2));
            end
        end
        function obj = redisc(obj, meta, varargin)
            
            switch meta.Name 
                case 'xPanels'
                    
                    xnew = linspace(0, obj.x(end,1), obj.xPanels);
                    
                case 'yPanels'
                    [rows, old] = size(obj.y);
                    new = obj.yPanels;
                    
                    % Point based therefore -1
                    old = (0:old-1)/(old - 1);
                    % Panel based
                    new = (0:new)/new;
                    
                    
                    for i = rows:-1:1
                        
                        xnew(i,:) = interp1(old, obj.x(i,:), new, 'pchip');
                        ynew(i,:) = interp1(old, obj.y(i,:), new, 'pchip');
                        znew(i,:) = interp1(old, obj.z(i,:), new, 'pchip');
                    end
            end
            
            obj.x = xnew;
            obj.y = ynew;
            obj.z = znew;
        end
        function obj = translate(obj, x, y, z)
            
            if nargin >= 2 && ~isempty(x)
                
                obj.x = obj.x + x;
            end
            if nargin >= 3 && ~isempty(y)
                
                obj.y = obj.y + y;
            end
            if nargin >= 4 && ~isempty(z)
                
                obj.z = obj.z + z;
            end
        end
        function obj = rotate(obj, rot)
            % Assume to be rotation around y-axis
            x = obj.x;
            z = obj.z;
            
            minx = min(x(:));
            obj.x = ((x - minx) * cos(rot) - z * sin(rot)) + minx;
            obj.z = z * cos(rot) + x * sin(rot);
        end
        function set.length(obj, val)
             
            obj.length = obj.radius * val;
        end
    end
end
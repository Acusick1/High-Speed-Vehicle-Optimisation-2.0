classdef Body < Geometry
    %AFTBODY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        
        height
        width
        area
        volume
    end
    
    properties
        conical = true
        %% TODO: Fixes for coarse bodies in wingbody
        % exits with ~any(nan) but interferes with first/last panels and errors
        yDisc = 0.5:0.025:1
        xDisc = 0.5*(1-cos(((0:80)*pi)./80))'
    end
    
    methods
        function a = get.width(obj)
            
            a = max(obj.y(:)) * 2;
        end
        function a = get.height(obj)
            
            zmid = obj.z(:, [end 1])';
            a = max(diff(zmid));
        end
        function a = get.area(obj)

            a = trapz(obj.x(:,1), max(obj.y, [], 2)); 
        end
        function a = get.volume(obj)
            
            % if isempty(obj.quad_data), obj.get_data; end
            % data = obj.quad_data;
            % a = 1/3 * sum(sum(dotmat(data.centre, data.unit_norm) .* data.area));
            
            % Calculates full volume (not half-body)
            [~,a] = convhull([obj.x(:); obj.x(:)], [obj.y(:); -obj.y(:)], ...
                [obj.z(:); obj.z(:)]);
        end
    end
    
    methods (Static)
        
        function [y, z] = combine_upper_lower(zu, zl, yu, yl)
            
            if nargin < 5 || isempty(yl), yl = yu; end
            
            y = [yu, fliplr(yl(:, 1:end-1))];
            z = [zu, fliplr(zl(:, 1:end-1))];
        end
        function [x, y, z] = close_body(x, y, z)
            
            [~, cols] = size(z);
            
            if ~all(z(1,:) == z(1))
                
                point = mean(z(1,:));
                z = [repmat(point, 1, cols); z];
                y = [zeros(1, cols); y];
                % Small offsets for unique interpolation purposes
                x = [x(1,:) - 1e-6; x];
            end
            
            if ~all(z(end,:) == z(end))
                
                point = mean(z(end,:));
                z = [z; repmat(point, 1, cols)];
                y = [y; zeros(1, cols)];
                x = [x; x(end,:) + 1e-6];
            end
        end
    end
end
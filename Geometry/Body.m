classdef Body < Geometry
    %AFTBODY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        conical = true;
    end
    
    methods
        
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
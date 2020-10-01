classdef Geometry < handle & matlab.mixin.Copyable 
    
    properties (GetObservable)
        
        points
        x
        y
        z
        sz
        tri_data
        quad_data
        update = true
    end
    
    methods
        function obj = Geometry()
            
            addlistener(obj,'x','PreGet',@obj.generate);
            addlistener(obj,'y','PreGet',@obj.generate);
            addlistener(obj,'z','PreGet',@obj.generate);
        end
        function obj = initialise(obj, varargin)
            
            obj.update = true;
        end
        function set.z(obj, val)
            
            obj.z = val;
            obj.sz = size(val);
        end
        function p = get.points(obj)
            
            p(:,:,1) = obj.x;
            p(:,:,2) = obj.y;
            p(:,:,3) = obj.z;
        end
        function obj = generate(obj, inputMetaProp, varargin)
            
            if isempty(obj.(inputMetaProp.Name)) || obj.update
                
                obj.update = false;
                obj.dogenerate();
%                 obj.quad_data = [];
%                 obj.tri_data = [];
            end
        end
        function obj = get_data(obj)
            
            quad.id = obj.quad_id(obj.sz);
            tri.id = obj.tri_id(quad.id);
            [quad.centre, quad.norm, quad.mag] = obj.quad_norm();
            unorm = quad.norm./quad.mag;
            % Collapsed panels have zero norm & mag therefore NaN
            unorm(isnan(unorm)) = 0;
            quad.unit_norm = unorm;
            
            tri.area = obj.tri_area(tri.id);
            quad.area = obj.quad_area(tri.area);
            
            obj.quad_data = quad;
            obj.tri_data = tri;
        end
        function area = tri_area(obj, id)
            
            X = obj.x(id);
            Y = obj.y(id);
            Z = obj.z(id);
            
            Dxy = X(:,1).*(Y(:,2) - Y(:,3)) - X(:,2).*(Y(:,1) - Y(:,3))...
                + X(:,3).*(Y(:,1) - Y(:,2));
            Dyz = Y(:,1).*(Z(:,2) - Z(:,3)) - Y(:,2).*(Z(:,1) - Z(:,3))...
                + Y(:,3).*(Z(:,1) - Z(:,2));
            Dzx = Z(:,1).*(X(:,2) - X(:,3)) - Z(:,2).*(X(:,1) - X(:,3))...
                + Z(:,3).*(X(:,1) - X(:,2));
            
            area = 0.5 * (Dxy.^2 + Dyz.^2 + Dzx.^2).^0.5;
        end
        function area = quad_area(obj, tri)
            
            if nargin < 2, tri = obj.tri_area(); end
            
            area = tri(1:2:end) + tri(2:2:end);
            area = reshape(area, obj.sz-1);
        end
        function [centre, norm, mag] = tri_norm(obj)
            
            id = obj.tri_id;
            
            p(:,:,1) = obj.x(id);
            p(:,:,2) = obj.y(id);
            p(:,:,3) = obj.z(id);
            
            centre = mean(p, 2);
            
            vec1 = p(:,2,:) - p(:,1,:);
            vec2 = p(:,3,:) - p(:,1,:);
            
            norm = zeros(size(vec1));
            
            % Ensure outward facing
            norm(1:2:end,:,:) = crossmat(vec1(1:2:end,:,:), vec2(1:2:end,:,:));
            norm(2:2:end,:,:) = crossmat(vec2(2:2:end,:,:), vec1(2:2:end,:,:));
            mag = magmat(norm);
        end
        function [centre, norm, mag] = quad_norm(obj)
            
            p = obj.points;
            
            a = p(1:end-1,1:end-1,:);
            b = p(2:end,1:end-1,:);
            c = p(2:end,2:end,:);
            d = p(1:end-1,2:end,:);
            
            centre = (a + b + c + d)/4;
            
            vec1 = c - a;
            vec2 = b - d;
            
            norm = crossmat(vec2, vec1);
            mag = magmat(norm);
        end
        function plot(obj, quad)
            
            if nargin < 2 || isempty(quad) || quad
            
                data = obj.quad_data;
            end
            
            cx = data.centre(:,:,1);
            cy = data.centre(:,:,2);
            cz = data.centre(:,:,3);
            
            nx = data.norm(:,:,1);
            ny = data.norm(:,:,2);
            nz = data.norm(:,:,3);
            
            figure
            hold on
            plot3(obj.x, obj.y, obj.z, 'k')
            plot3(obj.x', obj.y', obj.z', 'k')
            plot3(cx, cy, cz, 'k.')
            plot3(nx + cx, ny + cy, nz + cz, 'r.')
            axis equal
            hold off
        end
    end
    methods (Static)
        function id = quad_id(sz1, sz2)
            
            n = sz1(1);
            
            if nargin < 2 || isempty(sz2)
                
                m = sz1(2);
            else
                m = sz2;
            end
            
            mn = n * m;
            
            a = (1:mn - n)';
            
            delete = n:n:mn - n;
            a(delete,:) = [];
            
            id = [a, a+1, n+a, n+a+1];
        end
        function id = tri_id(qid)
            
            % if nargin < 1, qid = Geometry.quad_id; end
            
            [dim,~] = size(qid);
            
            id = zeros(dim*2,3);
            id(1:2:end,:) = qid(:,1:3);
            id(2:2:end,:) = qid(:,2:4);
        end
        function [init_obj, a] = define(init_obj, a)
            
            for i = 1:length(a)
                
                field = a(i).name;
                a(i).val = (a(i).min + a(i).max)/2;
                init_obj.(field) = a(i).val;
            end
        end
        function plotter(varargin)
            
            figure
            hold on
            for i = 1:length(varargin)
                
                if ~isempty(varargin{i})
                    
                    obj = varargin{i};
                    plot3(obj.x, obj.y, obj.z, 'k')
                    plot3(obj.x', obj.y', obj.z', 'k')
                end
            end
            view(0, 0)
            axis equal
            hold off
        end
    end
end
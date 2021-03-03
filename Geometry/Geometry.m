classdef Geometry% < handle & matlab.mixin.Copyable 
    
    properties %(SetObservable, GetObservable)
        
        points
        x
        y
        z
        centre
        sz
        tri_data
        quad_data
        update = true
        vio
    end
    
    properties (Dependent)
        
        nose_rad
        nose_id
    end
    
    methods
        function obj = Geometry()
            %% TODO: Quad_data not auto-updating
%             el(1) = addlistener(obj,'x','PreGet',@obj.generate);
%             el(2) = addlistener(obj,'y','PreGet',@obj.generate);
%             el(3) = addlistener(obj,'z','PreGet',@obj.generate);
%             
%             obj.el = el;
        end
        function obj = initialise(obj, varargin)
            
            obj.update = true;
        end
        function obj = set.z(obj, val)
            
            obj.z = val;
            obj = obj.get_data();
        end
        function a = get.sz(obj)
            
            a = size(obj.z);
        end
        function p = get.points(obj)
            
            p(:,:,1) = obj.x;
            p(:,:,2) = obj.y;
            p(:,:,3) = obj.z;
        end
        function obj = generate(obj, inputMetaProp, varargin)
            
            if obj.update || isempty(obj.(inputMetaProp.Name))
                
                obj.update = false;
                obj = obj.dogenerate();
                obj = obj.get_data();
                %% For testing
                % obj.update_plot();
            end
        end
        function obj = get_data(obj, varargin)
            
            qua.id = obj.quad_id(obj.sz);
            tri.id = obj.tri_id(qua.id);
            [qua.centre, qua.norm, qua.mag] = obj.quad_norm();
            unorm = qua.norm./qua.mag;
            % Collapsed panels have zero norm & mag therefore NaN
            unorm(isnan(unorm)) = 0;
            qua.unit_norm = unorm;
            
            tri.area = obj.tri_area(tri.id);
            qua.area = obj.quad_area(tri.area);
            
            obj.quad_data = qua;
            obj.tri_data = tri;
            
            %% TODO: Hack > Fix
            obj.centre = qua.centre;
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
            area(isnan(area)) = 0;
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
            
            %% Dealing with triangles
            
            t(:,:,1) = all(a == b, 3);
            t(:,:,2) = all(a == c, 3);
            t(:,:,3) = all(a == d, 3);
            t(:,:,4) = all(b == c, 3);
            t(:,:,5) = all(b == d, 3);
            t(:,:,6) = all(c == d, 3);
            
            % If two points and only two points are equal, panel is a
            % triangle
            t = sum(t, 3);
            istriangle = t == 1;
            
            if any(istriangle(:))
                
                dim = numel(istriangle);
                arr = (1:dim)';
                % Get [x, y, z] ids
                id = arr(istriangle) + (0:2) * dim;
                
                for i = size(id, 1):-1:1
                    
                    idi = id(i,:);
                    reduce = [a(idi); b(idi); c(idi); d(idi)];
                    keep = unique(reduce, 'rows', 'stable');
                    tri(i,:,:) = permute(keep, [3 2 1]);
                end
                
                % Replace tri panel data
                centre(id) = mean(tri, 3);
                vec1(id) = tri(:,:,3) - tri(:,:,1);
                vec2(id) = tri(:,:,2) - tri(:,:,1);
            end
            
            norm = crossmat(vec2, vec1);
            mag = magmat(norm);
        end
        function a = get.nose_id(obj)
           %% TODO: Hardcoding, correct nose location?
           % Also centre hack
            xNorm = obj.x - obj.x(1,:);
            xNorm = (xNorm(1:end-1, 1:end-1) + xNorm(1:end-1, 2:end) ...
                  + xNorm(2:end, 1:end-1) + xNorm(2:end, 2:end))/4;
            a = xNorm < 0.05;
        end    
        function a = get.nose_rad(obj)
            %% TODO: Hardcoding, correct nose location?
            % Translating to ensure nose start at (0,0)
            xInt = obj.x - obj.x(1,:);
            yInt = obj.y - obj.y(1,:);
            zInt = obj.z - obj.z(1,:);
            
            dim = size(xInt, 2);
            yNose = zeros(1, dim);
            zNose = zeros(1, dim);
            
            xDiff = diff(xInt, [], 1);
            
            %% Row removal (ensuring unique rows)
            row_con = all(xDiff == 0, 2);
            
            if any(row_con)
                
                xInt(row_con, :) = [];
                yInt(row_con, :) = [];
                zInt(row_con, :) = [];
                xDiff(row_con, :) = [];
            end
            
            %%
            % Matrix version of issorted, are columns sorted
            con = all(xDiff > 0, 1);
            
            if obj.conical
                
                % Assume constant x-distribution in radial dim
                if con(1)
                  
                    yNose = interp1(xInt(:,1), yInt, 0.05);
                    zNose = interp1(xInt(:,1), zInt, 0.05);
                end
            else
                for i = size(xInt, 2):-1:1
                    
                    % If xInt not sorted, interpolation will fail
                    if con(i)
                        
                        yNose(:,i) = interp1(xInt(:,i), yInt(:,i), 0.05);
                        zNose(:,i) = interp1(xInt(:,i), zInt(:,i), 0.05);
                    end
                end
            end
            
            rad = (yNose.^2 + zNose.^2).^0.5;
            
            % Fix for wings where chord almost perpendicular to x-axis
            % Should not be nan otherwise
            rad(isnan(rad)) = 1e-6;
            
            % Radius negative if x-norm is positive
            neg = obj.quad_data.norm(1,:,1) > 0;
            rad(neg) = -rad(neg);
            
            % Since forces etc calculated on panel centre, best to do same
            % with nose radius
            a = (rad(1:end-1) + rad(2:end))/2;
        end
        function splot = plot(obj, quad, surf_data)
            
            if nargin < 2 || isempty(quad) || quad
                
                data = obj.quad_data;
            end
            
            cx = data.centre(:,:,1);
            cy = data.centre(:,:,2);
            cz = data.centre(:,:,3);
            
            nx = data.norm(:,:,1);
            ny = data.norm(:,:,2);
            nz = data.norm(:,:,3);
            
            % Front and isometric
            views = [-90 0; -45, 25];
            % Top and bottom
            % views = [90 90; 90, -90];
            views = views(2,:);
            
            nviews = size(views, 1);
            figure(gcf)
            % clf
            for i = 1:nviews
                
                splot(i) = subplot(1, nviews, i);
                hold on
                if nargin < 3
                    plot3(obj.x, obj.y, obj.z, 'k')
                    plot3(obj.x', obj.y', obj.z', 'k')
                else
                    h(1) = surf(obj.x, obj.y, obj.z, surf_data);
                    h(2) = surf(obj.x, -obj.y, obj.z, surf_data);
                    set(h,'FaceLighting','flat','EdgeColor','none');
                    % colorbar
                end
                view(views(i,:))
                box on
                axis equal tight off
                hold off
            end
        end
        function update_plot(obj, varargin)
            
            if ishandle(1)
                
                figure(1)
                hold on
                clf
                obj.plot
            end
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
        function [a, var_names, nVar_parts] = init_variables(varargin)
            %% Initialise variables
            parts = varargin;
            
            a = [];
            for i = numel(parts):-1:1
                
                [~, b] = parts{i}(1).define;
                
                nVar_parts(i) = numel([b.val]);
                a = [b a];
            end
            
            nVar = sum(nVar_parts);
            var_names = repmat("", 1, nVar);
            
            i = 1;
            j = 1;
            while true
                
                dim = length(a(j).val);
                var_names((0:dim-1) + i) = a(j).name;
                i = i + dim;
                j = j + 1;
                if i > nVar, break; end
            end
        end
        function [init_obj, a] = define(init_obj, a)
            %% TODO: More generic? OptVariables?
            for i = 1:length(a)
                
                field = a(i).name;
                a(i).val = (a(i).min + a(i).max)/2;
                init_obj.(field) = a(i).val;
            end
        end
        function plotter(varargin)
            
            colour=[0.6 0.6 0.6];
            figure
            hold on
            for i = 1:length(varargin)
                
                if ~isempty(varargin{i})
                    
                    obj = varargin{i};
                    data = obj.quad_data;
                    cx = data.centre(:,:,1);
                    cy = data.centre(:,:,2);
                    cz = data.centre(:,:,3);
                    nx = data.unit_norm(:,:,1);
                    ny = data.unit_norm(:,:,2);
                    nz = data.unit_norm(:,:,3);
                    plot3(cx + nx, cy + ny, cz + nz, 'rx')
                    h = surf(obj.x,obj.y,obj.z,'LineWidth',0.05);
                    set(h,'FaceColor',colour,'FaceLighting','flat');%'EdgeColor','none');
                end
            end
           
            view(0, 0)
            xlabel('x, m', 'Interpreter', 'latex')
            ylabel('y, m', 'Interpreter', 'latex')
            zlabel('z, m', 'Interpreter', 'latex')
            plotFormat([], [], {'LineWidth', 0.5})
            axis equal tight
            hold off
        end
    end
end
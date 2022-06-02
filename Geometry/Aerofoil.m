classdef Aerofoil% < handle
    
    properties
        
        rescale = false
        xu = (0:0.01:1)'
        xl = (0:0.01:1)'
%         xu = 0.5*(1-cos(((0:0.01:1)*pi)))'
%         xl = 0.5*(1-cos(((0:0.01:1)*pi)))'
        conical = false
        checks
    end
    
    properties %(SetObservable)
        
        zu
        zl
    end
    
    properties (Dependent)
        
        data
        coords
        area
        thickness
        camber
        lead_edge
        trail_edge
    end
    
    properties (Constant)
        
        xLEr = 0.01
        % Minimum allowable thickness
        minThick = 1e-8;
        % Percentage of chord near LE/TE to remove from min thickness fixing
        edge = 0.05;
        maxTE = 0.05;
    end
    
    methods
        function obj = Aerofoil(varargin)
            
            if nargin >= 1
                if nargin == 1
                    
                    coords = varargin{:};
                    [~, front_id] = min(coords(:,1));
                    coords = coords - coords(front_id,:);
                    
                elseif nargin == 2
                    
                    if isboolen(varargin{2})
                        
                        [coords, obj.rescale] = deal(varargin{:});
                    else
                        [obj.zu, obj.zl] = deal(varargin{:});
                    end
                    
                elseif nargin >= 4
                    
                    [obj.zu, obj.zl, obj.xu, obj.xl] = ...
                        deal(varargin{1:4});
                    
                    if nargin == 5, obj.rescale = varargin{5}; end
                end
                
                if isempty(obj.zu)
                        
                    obj = obj.coordstoxz(coords);
                end
            end
            
            % addlistener(obj, 'zl', 'PostSet', @obj.check);
        end
        function data = get.data(obj)
            
            upNorm = Aerofoil.normal([obj.xu obj.zu]);
            loNorm = -Aerofoil.normal([obj.xl obj.zl]);
            
            norm(:,:,1) = [upNorm(:,1) loNorm(:,1)];
            norm(:,:,3) = [upNorm(:,2) loNorm(:,2)];
            
            mag_norm = magmat(norm);
            unit_norm = norm./mag_norm;
            points = reshape([obj.xu obj.xl obj.zu obj.zl], [], 2, 2);
            mag = magmat(diff(points, 1, 1));
            centre = (points(1:end-1,:,:) + points(2:end,:,:))/2;
            
            data.norm = norm;
            data.mag = mag_norm;
            data.unit_norm = unit_norm;
            data.area = mag;
            data.centre = centre;
        end
        function obj = coordstoxz(obj, coords)
            
            dim = ceil(length(coords)/2);
            % Split into upper and lower, ensure no repeating points
            upper = unique(coords(1:dim, :), 'rows');
            lower = unique(coords(dim:end, :), 'rows');
            
            %% Test > compare with aerofoilProperties hack
            upper([1 end], 1) = [0 1];
            lower([1 end], 1) = [0 1];
            
            %%
            if obj.rescale
                
                obj.xu = upper(:, 1);
                obj.xl = lower(:, 1);
                obj.zu = upper(:, 2);
                obj.zl = lower(:, 2);
            else
                if ~isequal(upper(:,1), obj.xu)
                    
                    obj.zu = obj.interp(upper);
                    % obj.zu = obj.interp(upper, [], [], 'pchip');
                end
                if ~isequal(lower(:,1), obj.xl)
                    
                    obj.zl = obj.interp(lower);
                    % obj.zl = obj.interp(lower, [], [], 'pchip');
                end
            end
        end
        function obj = redist(obj, x)
            
            x = x(:);
            
            zup = obj.interp(obj.xu, obj.zu, x, 'pchip');
            zlo = obj.interp(obj.xl, obj.zl, x, 'pchip');
            obj.xu = x;
            obj.xl = x;
            obj.zu = zup;
            obj.zl = zlo;
        end
        function obj = nondim(obj)
            %% TODO: nondim in x, translate in x and z
            minx = min(obj.coords(:,1));
            maxx = max(obj.coords(:,1));
            
            % Ensure non-dimensionalised
            obj.coords(:,1) = obj.coords(:,1) - minx;
            obj.coords = obj.coords/(maxx - minx);
        end
        function obj = dimensionalise(obj, c)
            
            obj.xu = obj.xu .* c;
            obj.xl = obj.xl .* c;
            obj.zu = obj.zu .* c;
            obj.zl = obj.zl .* c;
        end
        function obj = offset(obj, x, z)
            
            if nargin >= 2 && ~isempty(x)
                
                obj.xu = obj.xu + x;
                obj.xl = obj.xl + x;
            end
            if nargin >= 3 && ~isempty(z)
                
                obj.zu = obj.zu + z;
                obj.zl = obj.zl + z;
            end
        end
        function out = interp(obj, x, y, xq, varargin)
            
            if nargin < 3 || isempty(y)
                
                y = x(:,2);
                x = x(:,1);
            end
            
            if nargin < 4 || isempty(xq), xq = obj.xu; end
            if nargin < 5, varargin = {'linear', 'extrap'}; end
            
            try
                out = interp1(x, y, xq, varargin{:});
            catch
                [~, unique_id] = unique(x);
                x = x(unique_id,:);
                y = y(unique_id,:);
                
                out = interp1(x, y, xq, varargin{:});
            end
        end
        function obj = closeTrailEdge(obj)
            
            mid = (obj.zu(end) + obj.zl(end))/2;
            obj.zu(end) = mid;
            obj.zl(end) = mid;
        end
        function a = get.coords(obj)
            
            % 2:end avoids two leading edge points
            a = [flipud([obj.xu, obj.zu]);
                obj.xl(2:end), obj.zl(2:end)];
        end
        function a = check(obj, varargin)
            
            % varargin only here for event args or id
            
            [du_f, du_f2, ~, uMaxLoc] = obj.check_curve("upper");
            % [a(3), ~, ~, a(4)] = obj.check_curve("lower");
            [dt_f, ~, dt_r, tMaxLoc, tMax] = obj.check_curve("thickness");
            [~, tMin] = obj.check_thickness();
            LE = obj.lead_edge;
            LE = LE.radius;
            %a(12) = obj.thickness(end);
            
            a = [du_f, uMaxLoc, dt_f, dt_r, tMaxLoc, tMax, tMin, LE(1), ...
                LE(2), du_f2];
            
            if ~isempty(varargin) && isnumeric(varargin{1})
                
                id = varargin{1}; 
                a = a(id);
            end
        end
        function [d_f, d_f2, d_r, cLoc, Max] = check_curve(obj, which)
            
            switch which
                case "upper"
                    
                    x = obj.xu;
                    z = obj.zu;
                    
                case "lower"
                    
                    x = obj.xl;
                    z = -obj.zl;
                    
                case "thickness"
                    
                    x = obj.camber(:, 1);
                    z = obj.thickness;
            end
            
            % Calculates maximum non-dimensionalised aerofoil thickness
            [Max, ID] = max(z);
            % Location of maximum thickness
            cLoc = x(ID);
            
            % Gradients
            grad = diff(z)./diff(x);
            
            % Forward difference, diff returns vector sized length(x) - 1
            if ID == 1
                
                ID = 2;
            
            elseif ID == length(x)
                
                ID = ID - 1;
            end
            
            grad_f = grad(1:ID - 1);
            grad_f2 = diff(grad_f)./((diff(x(1:ID-1)) + diff(x(2:ID)))/2);
            grad_r = grad(ID:length(grad));
            
            d_f = min(grad_f);
            d_r = max(grad_r);
            
            if isempty(grad_f2)
                
                d_f2 = 0;
            else
                d_f2 = max(grad_f2);
            end
            
            if which == "lower", d_f = -d_f; d_r = -d_r; end
        end
        function [obj, tMin] = check_thickness(obj)
            
            c = obj.camber;
            t = obj.thickness;
            
            % Alter aerofoil if any thickness is < 0 or less than minimum thickness
            % outside of leading and trailing edges
            inner = c(:,1) > 0 & c(:,1) < 1;
            tFix = (t < obj.minThick & inner) | t < 0;

            obj = obj.fix_thickness(tFix);

            tMin = min(t(inner));
        end
        function obj = fix_thickness(obj, con)
            
            c = obj.camber;
            n = obj.normal(c);
            
            % Resetting thickness to min value perpendicular to camber line
            t = obj.minThick/2;
            
            cc = (c(1:end-1,:) + c(2:end,:))/2;
            cc = c;
            
            xu_int = cc(:,1) + t .* n(:,1);
            zu_int = cc(:,2) + t .* n(:,2);
            xl_int = cc(:,1) - t .* n(:,1);
            zl_int = cc(:,2) - t .* n(:,2);
            
            % Thickness applied in direction of camber, thus need to
            % interpolate back to orginal xu, xl locations
            obj.zu(con) = interp1(xu_int, zu_int, obj.xu(con), 'linear', 'extrap');
            obj.zl(con) = interp1(xl_int, zl_int, obj.xl(con), 'linear', 'extrap');
        end
        function a = get.area(obj)
            
            a = trapz((obj.xu + obj.xl)/2, obj.zu - obj.zl);
        end
        function a = get.thickness(obj)
            
            a = obj.zu - obj.zl;
        end
        function a = get.camber(obj)
            
            a(:,1) = (obj.xu + obj.xl)/2;
            a(:,2) = (obj.zu + obj.zl)/2;
        end
        function a = get.lead_edge(obj)
            
            % Derive leading edge radius, defined at 1% of chord
            zu001 = interp1(obj.xu, obj.zu, obj.xLEr, 'pchip');
            zl001 = interp1(obj.xl, obj.zl, obj.xLEr, 'pchip');
            
            a.radius = [zu001; -zl001];
            a.thickness = (zu001 - zl001)/2;
        end
        function a = get.trail_edge(obj)
            
            a = obj.zu(end,:) - obj.zl(end,:);
        end
        function plot(obj, figNum)
            
            if nargin >= 2 && isnumeric(figNum)
            
                f = figure(figNum);
            else
                f = figure;
            end
            
            ax = f.Children;
            
            if isempty(ax)
                
                nLines = 0;
            else
                nLines = numel(f.Children.Children);
            end
            
            lineSpec = {'k-', 'k--', 'k-.'};
            
            hold on
            plot(obj.coords(:,1), obj.coords(:,2), lineSpec{nLines+1})
            axis equal tight
            ylim([-0.2, 0.2])
            grid on
            xlabel('x/c', 'Interpreter', 'latex')
            ylabel('z/c', 'Interpreter', 'latex')
            f.Position = [50 50 850 420];
            
            if isa(obj, 'BezierFoil')
                
                plot(obj.control_points(:,1), obj.control_points(:,2), 'ko')
                plot(obj.control_points(:,3), obj.control_points(:,4), 'kx')
                
                if strcmp(obj.type, "tc")

                    figure(gcf)
                    hold on
                    plot(obj.camber(:,1), obj.thickness(:,1), 'k --') 
                    plot(obj.camber(:,1), obj.camber(:,2), 'k-.') 
                    hold off

                    legend('Aerofoil', 'Thickness Control Points', 'Camber Control Points', 'Thickness Curve', 'Camber Curve', 'interpreter', 'latex')
                else
                    legend('Aerofoil', 'Upper Control Points', 'Lower Control Points', 'interpreter', 'latex')
                end
            end
            
            plot_format()
            hold off
        end
    end
    
    methods (Static)
        function [aerofoil, file_names] = getaerofoil(names, path, id)
            
            if nargin < 1 || isempty(names)
                
                names = '.dat';
            end
            if nargin < 2 || isempty(path)
                
                path = fullfile(pwd, 'Geometry', 'DataFiles');
            end
            
            files = dir(path);
            file_names = convertCharsToStrings(...
                arrayfun(@(x) x.name, files, 'UniformOutput', false));
            
            file_names = file_names(contains(file_names, names, 'IgnoreCase', true));
            
            if nargin >= 3 && ~isempty(id)
                
                file_names = file_names(id);
            end
            
            for i = numel(file_names):-1:1
                
                try
                    str = fullfile(path, file_names(i));
                catch
                    str = join([path, file_names(i)], '/');
                end
                
                fid = fopen(str,'r');
                data = textscan(fid, '%f%f', 'HeaderLines', 1, 'Collect', 1);
                aerofoil(i,:) = Aerofoil(data{1});
                fclose(fid);
            end
        end
        function a = normal(c)
            
            [nx, nz] = normal2D(c(:,1), c(:,2));
            a = [nx, nz];
            %% TODO: Hack
            %a(end+1,:) = a(end,:);
        end
        function a = violation(id)
            
            val_min = [0, 0.15, -1e-6, nan, 0.15, 0.04, 0, 0.005, 0.001, nan];
            val_max = [nan, 0.65, nan, 1e-6, 0.65, 0.2, nan, nan, nan, 1e-3];
            norm = [0.5, nan, -1, -0.5, nan, nan, 0.01, 0.05, 0.05, 1e3];
            name = ["du_f", "cuMax", "dt_f", "dt_r", "ctMax", "tMax", "tMin", "rle", "rle", "du_f2"];

            if nargin < 1 || isempty(id), id = 1:numel(val_min); end
            
            a = Violation(val_min(id), val_max(id), norm(id), name(id));
        end
        function base = baseline(base_name)
            
            if nargin < 1
                
                base_name = "NACA66-206";
            end
            
            base = Aerofoil.getaerofoil(base_name);
            base.check;
        end
        function plotter(varargin)
            
            k = 1;
            
            for i = 1:numel(varargin)
                
                coords = [varargin{i}.coords];
                
                for j = 1:numel(varargin{i})
                    
                    figure
                    hold on
                    title(sprintf('Aerofoil%i', k))
                    plot(coords(:,j*2-1), coords(:,j*2), 'k');
                    axis equal
                    xlabel('x')
                    ylabel('z')
                    latex_figure()
                    hold off
                    k = k + 1;
                end
            end
        end
    end
end
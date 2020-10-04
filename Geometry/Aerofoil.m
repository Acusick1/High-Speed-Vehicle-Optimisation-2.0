classdef Aerofoil < handle
    
    properties
        
        rescale = false
        xu = (0:0.01:1)'
        xl = (0:0.01:1)'
        conical = false
        checks
    end
    
    properties (SetObservable)
        
        zu
        zl
    end
    
    properties (Dependent)
        
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
        tMin = 1e-3;
        % Percentage of chord near LE/TE to remove from min thickness fixing
        edge = 0.05;
    end
    
    methods
        function obj = Aerofoil(varargin)
            
            if nargin >= 1
                if nargin == 1
                    
                    coords = varargin{:};
                    
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
            
            addlistener(obj, 'zl', 'PostSet', @obj.check);
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
                end
                if ~isequal(lower(:,1), obj.xl)
                    
                    obj.zl = obj.interp(lower);
                end
            end
        end
        function obj = redist(obj, x)
            
            obj.zu = obj.interp(obj.xu, obj.zu, x);
            obj.zl = obj.interp(obj.xl, obj.zl, x);
            obj.xu = x;
            obj.xl = x;
        end
        function obj = nondim(obj)
            
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
        function out = interp(obj, x, y, xq)
            
            if nargin < 3 || isempty(y)
                
                y = x(:,2);
                x = x(:,1);
            end
            
            if nargin < 4, xq = obj.xu; end
            
            try
                out = interp1(x, y, xq, 'linear', 'extrap');
            catch
                [~, unique_id] = unique(x);
                x = x(unique_id,:);
                y = y(unique_id,:);
                
                out = interp1(x, y, xq, 'linear', 'extrap');
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
        function obj = check(obj, varargin)
            
            % varargin only here for event args
            [a(1), ~, a(2)] = obj.check_curve("upper");
            [a(3), ~, a(4)] = obj.check_curve("lower");
            [a(5), a(6), a(7), a(8)] = obj.check_curve("thickness");
            [obj, a(9)] = obj.check_thickness();
            LE = obj.lead_edge;
            a([10 11]) = LE.radius;
            obj.checks = a;
        end
        function [d_f, d_r, cLoc, Max] = check_curve(obj, which)
            
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
            grad_r = grad(ID:length(grad));
            
            if any(grad_f < 0)
                
                d_f = sum(grad_f(grad_f < 0));
            else
                d_f = min(grad_f);
            end
            
            if any(grad_r > 0)
                
                d_r = sum(grad_r(grad_r > 0));
            else
                d_r = max(grad_r);
            end
            
            if which == "lower", d_f = -d_f; d_r = -d_r; end
        end
        function [obj, out] = check_thickness(obj)
            
            c = obj.camber;
            t = obj.thickness;
            
            % Alter aerofoil if any thickness is < 0 or less than minimum thickness
            % outside of leading and trailing edges
            inner = c(:,1) > 0 & c(:,1) < 1;
            tFix = (t < obj.tMin & inner) | t < 0;

            obj = obj.fix_thickness(tFix);
          
            % If any negative values exist, provide minimum for penalty function
            % If no negative value exists, smallest thickness value outside of edges
            % will be returned. In this case, no penalty will be applied, but it
            % provides some continuity in the constraint function as opposed to always
            % returning a min value of zero (ie. at LE)
            if any(tFix)
                
                out = sum(t(tFix));
            else
                out = min(t);
            end
        end
        function obj = fix_thickness(obj, con)
            
            c = obj.camber;
            n = obj.normal(c);
            
            % Resetting thickness to min value perpendicular to camber line
            t = obj.tMin/2;
            
            cc = (c(1:end-1,:) + c(2:end,:))/2;
            
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
    end
    
    methods (Static)
        function [aerofoil, file_names] = getaerofoil(names)
            
            preloaded_dir = fullfile(pwd, 'Geometry', 'DataFiles');
            
            if nargin < 1
                
                names = 'dat';
            end
            
            files = dir(preloaded_dir);
            file_names = convertCharsToStrings(...
                arrayfun(@(x) x.name, files, 'UniformOutput', false));
            
            file_names = file_names(contains(file_names, names, 'IgnoreCase', true));
            
            for i = numel(file_names):-1:1
                
                try
                    str = fullfile(preloaded_dir, file_names(i));
                catch
                    str = join([preloaded_dir, file_names(i)], '/');
                end
                
                fid = fopen(str,'r');
                data = textscan(fid, '%f%f', 'HeaderLines', 1, 'Collect', 1);
                aerofoil = Aerofoil(data{1});
                fclose(fid);
            end
        end
        function a = normal(c)
            
            [nx, nz] = normal2D(c(:,1), c(:,2));
            a = [nx, nz];
        end
        function a = violation()
            
            a = Violation([0, 0.15, nan, 0.1, 0, 0.15, 0.05, 0, 0.005, 0.005], ...
                [nan, 0.65, 0, 0.65, nan, 0.65, 0.15, nan, nan, nan], [], ...
                [0.5, nan, -0.5, nan, nan, nan, nan, 0.01, 0.05, 0.05], ...
                ["du_f", "cuMax", "dl_f", "clMin", "dt_f", "ctMax", "tMax", "tMin", "rle", "rle"]);
        end
        function base = baseline(base_name)
            
            if nargin < 1
                
                base_name = "NACA66-206";
            end
            
            base = Aerofoil.getaerofoil(base_name);
            base.check;
            
%             a = violation([0, 0.15, nan, 0.1, 0, 0.15, 0.05, 0, 0.005, 0.005], ...
%                 [nan, 0.65, 0, 0.65, nan, 0.65, 0.15, nan, nan, nan], [], ...
%                 [0.5, nan, -0.5, nan, nan, nan, nan, 0.01, 0.05, 0.05], ...
%                 ["du_f", "cuMax", "dl_f", "clMin", "dt_f", "ctMax", "tMax", "tMin", "rle", "rle"]);
        end
        function plotter(varargin)
            
            k = 1;
            
            for i = 1:numel(varargin)
                
                coords = [varargin{i}.coords];
                
                for j = 1:numel(varargin{i})
                    
                    figure
                    hold on
                    title(sprintf('Aerofoil%i', k))
                    plot(coords(:,j*2-1), coords(:,j*2));
                    axis equal
                    k = k + 1;
                end
            end
        end
    end
end
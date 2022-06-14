classdef Aerofoil
    
    properties
        
        % Upper and lower aerofoil surface x and z coordinates
        xu = (0:0.01:1)'
        xl = (0:0.01:1)'
        % xu = 0.5*(1-cos(((0:0.01:1)*pi)))'
        % xl = 0.5*(1-cos(((0:0.01:1)*pi)))'
        zu
        zl
        % Interpolate input aerofoil coordinates to x values above
        rescale = false 
        conical = false         % 2D shape
        checks                  % Geometric violation parameters
        data
    end
    
    properties (Dependent)
        
        % Wrapped aerofoil coordinates from trailing edge to trailing edge
        coords
        area                    % Area of each 2D "panel"
        thickness               % Thickness of aerofoil at every x
        camber                  % Camber line of aerofoil
        lead_edge               
        trail_edge
    end
    
    properties (Constant)
        
        xLEr = 0.01             % Leading edge radius defined at 1% chord
        min_thickness = 1e-8;   % Minimum allowable thickness
    end
    
    methods
        function obj = Aerofoil(varargin)
            %AEROFOIL initialisation function
            %   Input options:
            %       Aerofoil(coords)
            %       Aerofoil(coords, rescale)
            %       Aerofoil(zu, zl)
            %       Aerofoil(zu, zl, xu, xl)
            %       Aerofoil(zu, zl, xu, xl, rescale)
            
            if nargin > 0
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
                        
                    obj = obj.coords2xz(coords);
                end
                
                obj = obj.get_data();
            end
        end
        function [obj, d] = get_data(obj)
            %% TODO: Should this be get function or just carried out during initialisation?
            %GETDATA calculates important aerofoil attributes from aerofoil
            %coordinates
            
            % Calculate outward facing normals on upper and lower surfaces
            [nx_upper, nz_upper] = normal2D(obj.xu, obj.zu);
            [nx_lower, nz_lower] = normal2D(obj.xl, obj.zl);
            
            normal(:,:,1) = [nx_upper nx_lower];
            normal(:,:,3) = [nz_upper nz_lower];
            
            % Normalise normals to unit normals
            mag_norm = magmat(normal);
            unit_norm = normal./mag_norm;
            
            points(:,:,1) = [obj.xu obj.xl];
            points(:,:,2) = [obj.zu obj.zl];
            
            % "area" of line segment = length of line segment
            line_area = magmat(diff(points, 1, 1));
            
            % Get centre point of each line segment
            centre = (points(1:end-1,:,:) + points(2:end,:,:))/2;
            
            d.norm = normal;
            d.mag = mag_norm;
            d.unit_norm = unit_norm;
            d.area = line_area;
            d.centre = centre;
            
            obj.data = d;
        end
        function obj = coords2xz(obj, coords)
            %COORDS2XZ sets upper and lower leading to trailing edge
            %coordinates from input wrapped coordinates
            %   Inputs:
            %   coords - wrapped aerofoil coordinates with start/end points 
            %      at trailing edge
            
            nPoints = ceil(length(coords)/2);
            % Split into upper and lower, ensure no repeating points
            upper = unique(coords(1:nPoints, :), 'rows');
            lower = unique(coords(nPoints:end, :), 'rows');
            
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
            %REDIST redistributes current aerofoil by input x values
            %   Inputs:
            %   x - non-dimensionalised array of increasing x values from 0
            %       to 1
            
            x = x(:);
            
            zup = obj.interp(obj.xu, obj.zu, x, 'pchip');
            zlo = obj.interp(obj.xl, obj.zl, x, 'pchip');
            obj.xu = x;
            obj.xl = x;
            obj.zu = zup;
            obj.zl = zlo;
        end
        function obj = nondim(obj)
            %NONDIM translate and rescale aerofoil coordinates to
            %non-dimensionalised coordinates with a chord 0 to 1
            
            % Translate x and z coordinates by leading edge value
            [~, lead_edge_id] = min(obj.coords(:,1));
            obj.coords = obj.coords - obj.coords(lead_edge_id,:);
            
            % Leading edge now [0, 0], so divide by maximum chord value
            % (trailing edge x) in both x and z to non-dimensionalise
            [~, trail_edge_id] = max(obj.coords(:,1));
            obj.coords = obj.coords/obj.coords(trail_edge_id, 1);
        end
        function obj = dimensionalise(obj, chord)
            %DIMENSIONALISE transforms non-dimensionalised aerofoil 
            %coordinates by input scalar chord
            
            obj.xu = obj.xu * chord;
            obj.xl = obj.xl * chord;
            obj.zu = obj.zu * chord;
            obj.zl = obj.zl * chord;
        end
        function obj = offset(obj, x, z)
            %OFFSET traslates non-dimensional aerofoil chord by input 
            %x and z scalar values
            
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
            %INTERP aerofoil specific interpolation function
            %   wrapper function of interp1 which requires unique points,
            %   not guaranteed in aerofoil generation
            %   
            %   Inputs: equivalent to interp1 function
            
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
        function obj = close_trail_edge(obj)
            %CLOSE_TRAIL_EDGE transforms finite thickness trailing edge to
            %point
            mid = (obj.zu(end) + obj.zl(end))/2;
            obj.zu(end) = mid;
            obj.zl(end) = mid;
        end
        function params = check(obj, varargin)
            %% TODO: does this need varargin or only for handle implementation?
            % varargin only here for event args or id
            
            [du_f, du_f2, ~, uMaxLoc] = obj.check_curve("upper");
            % [a(3), ~, ~, a(4)] = obj.check_curve("lower");
            [dt_f, ~, dt_r, tMaxLoc, tMax] = obj.check_curve("thickness");
            [~, tMin] = obj.constrain_thickness();
            LE = obj.lead_edge;
            LE = LE.radius;
            %a(12) = obj.thickness(end);
            
            params = [du_f, uMaxLoc, dt_f, dt_r, tMaxLoc, tMax, tMin, ...
                LE(1), LE(2), du_f2];
            
            if ~isempty(varargin) && isnumeric(varargin{1})
                
                id = varargin{1}; 
                params = params(id);
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
        function [obj, min_thick] = constrain_thickness(obj)
            %CONSTRAIN_THICKNESS ensures aerofoil thickness is greater than
            %min_thickness
            
            % Alter aerofoil if any thickness is < 0 or less than minimum 
            % thickness outside of leading and trailing edges
            inner = obj.camber(:,1) > 0 & obj.camber(:,1) < 1;
            
            fix_thickness = ...
                (obj.thickness < obj.min_thickness & inner) | ...
                obj.thickness < 0;
            
            if any(fix_thickness)
                obj = obj.fix_thickness(fix_thickness);
            end
            
            % Return minimum thickess of constrained aerofoil
            min_thick = min(obj.thickness(inner));
        end
        function obj = fix_thickness(obj, fix)
            %FIX_THICKNESS sets aerofoil thickness less than min_thickness
            %to min_thickness
            %   Inputs:
            %   fix - logical array representing points that are below 
            %       min_thickness 
            
            [nx, nz] = normal2D(obj.camber(:,1), obj.camber(:,2));
            % Taking central difference since normal outputs n-1 elements
            c = (obj.camber(1:end-1,:) + obj.camber(2:end,:))/2;
            
            % Resetting thickness to min value perpendicular to camber line
            t = obj.min_thickness/2;
            
            xu_int = c(:,1) + t .* nx;
            zu_int = c(:,2) + t .* nz;
            xl_int = c(:,1) - t .* nx;
            zl_int = c(:,2) - t .* nz;
            
            % Thickness applied in direction of camber, thus need to
            % interpolate back to orginal xu, xl locations
            obj.zu(fix) = interp1(xu_int, zu_int, obj.xu(fix), 'linear', 'extrap');
            obj.zl(fix) = interp1(xl_int, zl_int, obj.xl(fix), 'linear', 'extrap');
        end
        function coords = get.coords(obj)
            
            if ~isempty(obj.zu) && ~isempty(obj.zl)
                % 2:end avoids two leading edge points
                coords = [flipud([obj.xu, obj.zu]);
                          obj.xl(2:end), obj.zl(2:end)];
            end
        end
        function area = get.area(obj)
            
            area = trapz((obj.xu + obj.xl)/2, obj.zu - obj.zl);
        end
        function thickness = get.thickness(obj)
            
            thickness = obj.zu - obj.zl;
        end
        function camber = get.camber(obj)
            
            camber(:,1) = (obj.xu + obj.xl)/2;
            camber(:,2) = (obj.zu + obj.zl)/2;
        end
        function lead_edge = get.lead_edge(obj)
            
            % Derive leading edge radius, defined at 1% of chord
            zu001 = interp1(obj.xu, obj.zu, obj.xLEr, 'pchip');
            zl001 = interp1(obj.xl, obj.zl, obj.xLEr, 'pchip');
            
            lead_edge.radius = [zu001; -zl001];
            lead_edge.thickness = (zu001 - zl001)/2;
        end
        function trail_edge = get.trail_edge(obj)
            
            trail_edge = obj.zu(end,:) - obj.zl(end,:);
        end
        function plot(obj, fig_num)
            %% TODO: Still relevant? Only has BezierFoil option
            
            if nargin >= 2 && isnumeric(fig_num)
            
                f = figure(fig_num);
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
        function [aerofoil, file_names] = get_aerofoil(pattern, path)
            %GET_AEROFOIL loads aerofoil data to consistent format
            %   Creates Aerofoil objects from data files
            %
            %   Inputs:
            %   pattern - filename of aerofoil dat file, or pattern to load
            %       multiple aerofoils (default = all dat files), string or
            %       cell arrays of characters also allowed
            %   path - path to directory containing data files
            %       (default = "BasePath/Geometry/DataFiles")
            
            if nargin < 1 || isempty(pattern)
                
                pattern = '.dat';
            end
            if nargin < 2 || isempty(path)
                
                path = fullfile(get_base_path(), 'Components', 'Geometric', 'DataFiles');
            end
            
            % Get files within path, extract all names
            files = dir(path);
            file_names = convertCharsToStrings(...
                arrayfun(@(x) x.name, files, 'UniformOutput', false));
            
            % Narrow names down to those that match input patterns
            file_names = file_names(...
                contains(file_names, pattern, 'IgnoreCase', true));
            
            for i = numel(file_names):-1:1
                
                fid = fopen(fullfile(path, file_names(i)), 'r');
                data = textscan(...
                    fid, '%f%f', 'HeaderLines', 1, 'Collect', 1);
                
                % Create object from loaded data
                aerofoil(i,:) = Aerofoil(data{1});
                fclose(fid);
            end
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
            
            base = Aerofoil.get_aerofoil(base_name);
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
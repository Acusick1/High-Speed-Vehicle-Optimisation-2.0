classdef BezierFoil < Aerofoil & Bezier
    
    properties
        
        type = "direct"
        bez
        names
    end
    
    methods
        function [obj] = BezierFoil(control_points, type)
            
            addlistener(obj,'control_points','PostSet',@obj.dogenerate);
            
            if nargin > 0
                %% TODO: Addlistener in superclass? generate > dogenerate 
                obj.control_points = control_points;
                
                %%
                if nargin > 3 && ~isempty(type), obj.type = type; end
                
                switch obj.type
                
                    case "direct"

                        vars = ["cpxu" "cpzu" "cpxl" "cpzl"];

                    case "tc"

                        vars = ["cpxc" "cpzc" "cpxt" "cpzt"];
                end
                
                obj.set_names(vars);
            end
        end
        function obj = dogenerate(obj, varargin)
            
            % Make Bezier curves from control points
            %% TODO
            % combining these to be made as array of objects results in
            % objects being the same since class was changed from value to
            % handle
            new_bez = ...
                Bezier(obj.control_points(:,1:2), obj.control_points(:,3:4));
            %%
            
            switch obj.type
                case "tc"
                    % Bezier curves equal aerofoil thickness and camber lines
                    t = new_bez(1).curve;
                    c = new_bez(2).curve;
                    
                    t(:,2) = interpulinex(t, c(:,1));
                    n = Aerofoil.normal(c);
                    
                    xu = c(:,1) + t(:,2) .* n(:,1);
                    zu = c(:,2) + t(:,2) .* n(:,2);
                    xl = c(:,1) - t(:,2) .* n(:,1);
                    zl = c(:,2) - t(:,2) .* n(:,2);
                    
                case "direct"
                    
                    % Bezier curves represent upper and lower aerofoil surfaces
                    upper = new_bez(1).curve;
                    lower = new_bez(2).curve;
                    
                    xu = upper(:,1);
                    zu = upper(:,2);
                    xl = lower(:,1);
                    zl = lower(:,2);
            end
            
            obj.zu = obj.interp(xu, zu);
            obj.zl = obj.interp(xl, zl);
            obj.bez = new_bez;
        end
        function set_names(obj, vars)
            
            [n,~] = size(obj.control_points);
            
            leng = n * length(vars);
            a = strings(leng, 1);
            
            k = 1;
            for i = 1:length(vars)
                for j = n:-1:1
                    
                    a(k) = vars(i) + num2str(j);
                    k = k + 1;
                end
            end
            
            obj.names = a;
        end
    end
    
    methods (Static)
        function aerofoil = test(type)
            
            if nargin < 1, type = "direct"; end
            
            switch type
                case "direct"
                    
                    cp1(:,1) = [1, 0.75, 0.6, 0.4, 0.29, 0.24, 0];
                    cp1(:,2) = [0, 0.01, 0.09, -0.03, 0, 0.11, 0];
                    cp2(:,1) = [1, 0.75, 0.6, 0.4, 0.5, 0.05, 0];
                    cp2(:,2) = [0, 0, -0.05, -0.06, -0.06, -0.07, 0];
                    
                case "tc"
                    
                    cp1(:,1) = [1, 0.75, 0.6, 0.4, 0.29, 0.24, 0];
                    cp1(:,2) = [0, 0.01, 0.05, 0.1, 0.05, 0.1, 0];
                    cp2(:,1) = [1, 0.75, 0.6, 0.5, 0.3, 0.05, 0];
                    cp2(:,2) = [0, 0, 0.01, 0.02, 0.02, 0.01, 0];
            end
            
            b = [cp1 cp2];
            
            [aerofoil] = BezierFoil(b, type);
            aerofoil.plotter(aerofoil)
            aerofoil.bez(1).plot
            aerofoil.bez(2).plot
        end
        function [init, a] = define(type, n, crossover, zMean, zVar, edgeVar)
            
            % Generic properties
            if nargin < 1 || isempty(type), type = "direct"; end
            if nargin < 2 || isempty(n), n = 6; end
            if nargin < 3 || isempty(crossover), crossover = 0; end
            if nargin < 4 || isempty(zMean), zMean = 0.05; end
            if nargin < 5 || isempty(zVar), zVar = 0.15; end
            if nargin < 6 || isempty(edgeVar), edgeVar = 0; end
            
            % Front and end points fixed
            nFree = n - 2;
            
            x = 1:-1/nFree:0;
            xMin = [1, max(x - 1/nFree - crossover, 0)]';
            xMax = [min(x + 1/nFree + crossover, 1) ,0]';
            
            zlMin = [-edgeVar, repmat(-zMean - zVar, 1, nFree), 0]';
            zlMax = [edgeVar, repmat(-zMean + zVar, 1, nFree), 0]';
            
            zuMin = [-edgeVar, repmat(zMean - zVar, 1, nFree), 0]';
            zuMax = [edgeVar, repmat(zMean + zVar, 1, nFree), 0]';
            
            var_min = [xMin, zuMin, xMin, zlMin];
            var_max = [xMax, zuMax, xMax, zlMax];
            
            init = (var_min + var_max)/2;
            
            init = BezierFoil(init, type);
            names = repmat("control_points", size(var_min));
            
            if crossover, init.con_xcp = true; end
            
            % a = OptVariables(Min, Max, names, con, trans, opt);
            a.name = names(:)';
            a.min = var_min(:)';
            a.max = var_max(:)';
            a.val = (var_min(:) + var_max(:))'/2;
        end
    end
end
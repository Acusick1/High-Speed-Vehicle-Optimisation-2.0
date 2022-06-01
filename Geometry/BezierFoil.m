classdef BezierFoil < Aerofoil & Bezier
    
    properties
        
        type = "direct"
        names
    end
    
    methods
        function [obj] = BezierFoil(control_points, type)
            
            % addlistener(obj,'control_points','PostSet',@obj.dogenerate);
            
            if nargin > 0
                
                if nargin >= 2 && ~isempty(type), obj.type = type; end
                % Here for set method, relies on arg 2
                obj.control_points = control_points;
                
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
            
            switch obj.type
                case "tc"
                    % Bezier curves equal aerofoil thickness and camber lines
                    t = obj.curve(:,:,1);
                    c = obj.curve(:,:,2);
                    
                    t(:,2) = interp1(t(:,1), t(:,2), c(:,1));
                    n = Aerofoil.normal(c);
                    
                    xu = c(:,1) + t(:,2) .* n(:,1);
                    zu = c(:,2) + t(:,2) .* n(:,2);
                    xl = c(:,1) - t(:,2) .* n(:,1);
                    zl = c(:,2) - t(:,2) .* n(:,2);
                    
                case "direct"
                    
                    % Bezier curves represent upper and lower aerofoil surfaces
                    upper = obj.curve(:,:,1);
                    lower = obj.curve(:,:,2);
                    
                    xu = upper(:,1);
                    zu = upper(:,2);
                    xl = lower(:,1);
                    zl = lower(:,2);
            end
            
            obj.zu = obj.interp(xu, zu);
            obj.zl = obj.interp(xl, zl);
        end
        function obj = set_names(obj, vars)
            
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
                    
                    cp1(:,1) = [1, 0.75, 0.6, 0.4, 0.29, 0.14, 0];
                    cp1(:,2) = [0, 0.01, 0.09, 0.06, 0.02, 0.11, 0];
                    cp2(:,1) = [1, 0.75, 0.6, 0.4, 0.5, 0.05, 0];
                    cp2(:,2) = [0, 0, -0.05, -0.06, -0.06, -0.07, 0];
                    
                case "tc"
                    
                    cp1(:,1) = [1, 0.75, 0.6, 0.4, 0.25, 0.1, 0];
                    cp1(:,2) = [0, 0.01, 0.05, 0.08, -0.025, 0.06, 0];
                    cp2(:,1) = [1, 0.75, 0.6, 0.5, 0.3, 0.05, 0];
                    cp2(:,2) = [0, 0, 0.01, 0.02, 0.05, 0.03, 0];
            end
            
            b = [cp1 cp2];
            
            [aerofoil] = BezierFoil(b, type);
            aerofoil.plot(aerofoil)
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
            
            xlin = 1:-1/nFree:0;
            xhcos = 1-cos((xlin*(pi/2))/max(xlin));
            
            xMin = [1, max(xhcos - 1/nFree - crossover, 0)]';
            xMax = [min(xhcos + 1/nFree + crossover, 1) ,0]';
            
            if strcmpi(type, "direct")
                
                zlMin = [-edgeVar/2, repmat(-zMean - zVar, 1, nFree), 0]';
                zlMax = [edgeVar/2, repmat(-zMean + zVar, 1, nFree), 0]';
                zlMax = min(zlMax, 0);
                
                zuMin = [-edgeVar/2, repmat(zMean - zVar, 1, nFree), 0]';
                zuMax = [edgeVar/2, repmat(zMean + zVar, 1, nFree), 0]';
                zuMin = max(zuMin, 0);
                
                var_min = [xMin, zuMin, xMin, zlMin];
                var_max = [xMax, zuMax, xMax, zlMax];
                
            elseif strcmpi(type, "tc")
            
                ztMin = [0, repmat(zMean - zVar, 1, nFree), 0]';
                ztMax = [edgeVar, repmat(zMean + zVar, 1, nFree), 0]';
                
                zcMin = [-edgeVar/2, repmat(0 - zVar, 1, nFree), 0]';
                zcMax = [edgeVar/2, repmat(zMean + zVar, 1, nFree), 0]';
                
                var_min = [xMin, ztMin, xMin, zcMin];
                var_max = [xMax, ztMax, xMax, zcMax];
            end
            
            var_init = (var_min + var_max)/2;
            
            init = BezierFoil(var_init, type);
            names = repmat("control_points", size(var_min));
            
            if crossover, init.con_xcp = true; end
            if strcmpi(type, "direct"), init.con_zend = true; end
            
            % a = OptVariables(Min, Max, names, con, trans, opt);
            a.name = names(:)';
            a.min = var_min(:)';
            a.max = var_max(:)';
            a.val = (var_min(:) + var_max(:))'/2;
        end
    end
end
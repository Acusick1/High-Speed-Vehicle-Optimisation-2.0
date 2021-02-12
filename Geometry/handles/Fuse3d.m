classdef Fuse3d < Body
    
    properties (SetObservable)
        
        length
        w
        h
        nc
        nd
        Acu
        Acl
        Ad
        upper
        lower
    end
    properties
        
        shapefuns
        classfuns
    end
    
    methods
        function obj = Fuse3d(length, nc, nd, w, h, Acu, Acl, Ad)
            
            if nargin > 0
                
                obj.length = length;
                % Halved as defining half width and upper/lower 
                obj.w = w * length * 0.5;
                obj.h = h * length * 0.5;
                obj.nc = nc;
                obj.nd = nd;
                obj.Acu = Acu;
                obj.Acl = Acl;
                obj.Ad = Ad;
            end
            
            addlistener(obj,'length','PostSet',@obj.initialise);
            addlistener(obj,'w','PostSet',@obj.initialise);
            addlistener(obj,'h','PostSet',@obj.initialise);
            addlistener(obj,'nc','PostSet',@obj.initialise);
            addlistener(obj,'nd','PostSet',@obj.initialise);
            addlistener(obj,'Acu','PostSet',@obj.initialise);
            addlistener(obj,'Acl','PostSet',@obj.initialise);
            addlistener(obj,'Ad','PostSet',@obj.initialise);
        end
        function obj = dogenerate(obj)
            
            S = obj.shapefuns;
            C = obj.classfuns;
            
            y_nd = -S.d .* C.d .* (1 - (2 * obj.yDisc));
            zu_nd = obj.h(1) * C.cu .* S.cu .* (C.d .* S.d);
            zl_nd = obj.h(2) * C.cl .* S.cl .* (C.d .* S.d);
            
            % Max or max - min as min will always be 0
            y = y_nd./max(y_nd(:)) * obj.w * obj.length * 0.5;
            zu = zu_nd./max(zu_nd(:)) * obj.h(1) * obj.length * 0.5;
            zl = zl_nd./max(-zl_nd(:)) * obj.h(end) * obj.length * 0.5;
            
            [y, z] = obj.combine_upper_lower(zu, zl, y);
            
            x = obj.xDisc * obj.length;
            x = repmat(x, 1, size(y, 2));
            
            [obj.x, obj.y, obj.z] = obj.close_body(x, y, z);
            
            obj.upper = zu;
            obj.lower = zl;
            obj.get_data;
        end
        function set.Acu(obj, in)
            %% TODO: Yay or nae
            % Can only be slightly below previous
%             for i = 2:numel(in)
%                 
%                 if in(i) < in(i-1)*0.9
%                     
%                     in(i) = in(i-1)*0.9;
%                 end
%             end
            
            obj.Acu = in; 
        end
        function set.Acl(obj, in)
            
            % Can only be slightly above previous
%             for i = 2:numel(in)
%                 
%                 if in(i) > in(i-1) + 0.1
%                     
%                     in(i) = in(i-1) + 0.1;
%                 end
%             end
            
            obj.Acl = in; 
        end
        function a = check(obj)
            
            upper_diff = obj.upper(:,2:end) - obj.upper(:,1);
            lower_diff = -(obj.lower(:,2:end) - obj.lower(:,1));
            
            a = [max(upper_diff(:)), max(lower_diff(:)), ...
                mean(obj.nose_rad)/obj.length, obj.height/obj.width];
        end
        function a = get.shapefuns(obj)
            
            Au = obj.Acu;
            Al = obj.Acl;
            
            a.cu = obj.shapefun([Au, flip(Au)], obj.yDisc);
            a.cl = obj.shapefun([Al, flip(Al)], obj.yDisc);
            a.d = obj.shapefun(obj.Ad, obj.xDisc);
        end
        function a = get.classfuns(obj)
            
            a.cu = obj.classfun(obj.nc(1), obj.yDisc);
            a.cl = obj.classfun(obj.nc(2), obj.yDisc);
            a.d = obj.classfun(obj.nd, obj.xDisc);
        end
        function out = interp_shape(obj, x, y)
            
            funs = obj.shapefuns;
            
            if nargin < 1 || isempty(x)
                
                out(:,1) = funs.Sd;
            else
                out(:,1) = interp1(obj.xDisc, funs.Sd, x);
            end
            
            if nargin < 2 || isempty(y)
                
                out(:,2) = funs.Scu;
                out(:,3) = funs.Scl;
            else
                out(:,2) = interp1(obj.yDisc, funs.Scu, y);
                out(:,3) = interp1(obj.yDisc, funs.Scl, y);
            end
        end
        function out = interp_class(obj, x, y)
            
            funs = obj.classfuns;
            
            if nargin < 1 || isempty(x)
                
                out(:,1) = funs.Cd;
            else
                out(:,1) = interp1(obj.xDisc, funs.Cd, x);
            end
            
            if nargin < 2 || isempty(y)
                
                out(:,2) = funs.Scu;
                out(:,3) = funs.Scl;
            else
                out(:,2) = interp1(obj.yDisc, funs.Ccu, y);
                out(:,3) = interp1(obj.yDisc, funs.Ccl, y);
            end
        end
        function plotbody(obj)
            
            Geometry.plotter(obj)
            figure(gcf)
            hold on
            title(sprintf(...
                "NC = [%.2f, %.2f], ND = [%.2f, %.2f]", obj.nc, obj.nd));
            hold off
        end
        function plotcurves(obj)
            
            S = obj.shapefuns;
            C = obj.classfuns;
            
            y = -S.d .* C.d;
            zu = C.cu .* S.cu;
            zl = C.cl .* S.cl;
            
%             Geometry.plotter(obj)
            yspare = zeros(size(y));
            zspare = zeros(size(zu));
            figure(gcf)
            hold on
            plot3(obj.xDisc, -y, yspare)
            plot3(zspare, obj.yDisc-0.5, zu)
            plot3(zspare, obj.yDisc-0.5, zl)
            hold off
        end
    end
    methods (Static)
        function S = shapefun(A, disc)
            
            n = length(A) - 1;
            S = 0;
            
            for i = 0:n
                
                K = nchoosek(n, i);
                
                Si = K * disc.^i .* (1 - disc).^(n-i);
                S = S + A(i+1) * Si;
            end
        end
        function C = classfun(n, disc, scale)
            
            if length(n) == 1, n(2) = n(1); end
            if nargin < 3 || isempty(scale), scale = 1; end
            
            C = scale * disc.^n(1) .* (1 - disc).^n(2);
        end
        function [init_obj, a] = define()

            a(1) = struct('name', "length", 'min', 1, 'max', 20);
            % Stretch factors
            a(2) = struct('name', "w", 'min', 0.1, 'max', 1);
            a(3) = struct('name', "h", 'min', [0.1 0.1], 'max', [1 1]);
            % Soanwise upper/lower curvature: 0 = tall/box, 1 = flat/sharp
            a(4) = struct('name', "nc", 'min', [0.05 0.05], 'max', [1 1]);
            % Streamwise nose/tail curvature: 0 = flat, 1 = knife
            a(5) = struct('name', "nd", 'min', [0.25 1e-3], 'max', [1 1]);
            %% TODO: linear Ac1 > Ac2 rather than fuse height constraint
            a(6) = struct('name', "Acu", 'min', [0.05 0.05], 'max', [1 1]);
            a(7) = struct('name', "Acl", 'min', -[1 1], 'max', -[0.02 0.02]);
            % Have to be > 0 otherwise all(z) == 0
            a(8) = struct('name', "Ad", 'min', [1e-3 1e-3 1e-3 1e-3], 'max', [1 1 1 1]);            
            
            init_obj = Fuse3d();
            [init_obj, a] = Geometry.define(init_obj, a);
        end
        function obj = test()
            
            L = 3.5;
            c = [4, 5];
            d = 5;
            nc = [0.949, 0.949];
            nd = [0.542, 0.105];
            Acu = flip([0.462, 0.136]);
            Acl = -[0.044, 0.371];
            Ad = [0.258, 0.159, 0.170, 0.199];
            obj = Fuse3d(L, nc, nd, c, d, Acu, Acl, Ad);
            obj.dogenerate;
            obj.plot()
        end
    end
end


classdef Fuse3d < Body
    
    properties (SetObservable)
        
        length
        nc
        nd
        c
        d
        Acu
        Acl
        Ad
    end
    properties

        height
        width
        area
        volume
        shapefuns
        classfuns
    end
    properties (Constant)
        
        name = "body"
        %% TODO: Fixes for coarse bodies in wingbody
        % exits with ~any(nan) but interferes with first/last panels and errors
        xDisc = (0:0.01:1)'
        yDisc = 0.5:0.01:1
    end
    
    methods
        function obj = Fuse3d(length, nc, nd, c, d, Acu, Acl, Ad)
            
            if nargin > 0
                
                if nargin < 6 || isempty(c), c = [1, 1]; end
                if nargin < 7 || isempty(d), d = 1; end
                
                obj.length = length;
                obj.c = c;
                obj.d = d;
                obj.nc = nc;
                obj.nd = nd;
                obj.Acu = Acu;
                obj.Acl = Acl;
                obj.Ad = Ad;
            end
            
            addlistener(obj,'length','PostSet',@obj.initialise);
            addlistener(obj,'c','PostSet',@obj.initialise);
            addlistener(obj,'d','PostSet',@obj.initialise);
            addlistener(obj,'nc','PostSet',@obj.initialise);
            addlistener(obj,'nd','PostSet',@obj.initialise);
            addlistener(obj,'Acu','PostSet',@obj.initialise);
            addlistener(obj,'Acl','PostSet',@obj.initialise);
            addlistener(obj,'Ad','PostSet',@obj.initialise);
        end
        function obj = dogenerate(obj)
            
            S = obj.shapefuns;
            C = obj.classfuns;
            
            y = -S.d .* C.d .* (1 - (2 * obj.yDisc));
            zu = C.cu .* S.cu .* (C.d .* S.d);
            zl = C.cl .* S.cl .* (C.d .* S.d);
            
            [y, z] = obj.combine_upper_lower(zu, zl, y);
            
            x = obj.xDisc * obj.length;
            x = repmat(x, 1, size(y, 2));
            
            [obj.x, obj.y, obj.z] = obj.close_body(x, y, z);
            
            obj.height = max(zu(:,1) - zl(:,1));
            obj.width = max(y(:)) * 2;
            obj.area = trapz(x(:,1), max(y, [], 2));
            obj.get_data;
        end
        function a = get.volume(obj)
            
            data = obj.quad_data;
            a = 1/3 * sum(sum(dotmat(data.centre, data.unit_norm) .* data.area));
        end
        function a = get.height(obj)
            
            zmid = obj.z(:, [end 1])';
            a = max(diff(zmid));
        end
        function a = get.shapefuns(obj)
            
            Au = obj.Acu;
            Al = obj.Acl;
            
            a.cu = obj.shapefun([Au, flip(Au)], obj.yDisc);
            a.cl = obj.shapefun([Al, flip(Al)], obj.yDisc);
            a.d = obj.shapefun(obj.Ad, obj.xDisc);
        end
        function a = get.classfuns(obj)
            
            a.cu = obj.classfun(obj.nc(1), obj.yDisc, obj.c(1));
            a.cl = obj.classfun(obj.nc(2), obj.yDisc, obj.c(2));
            a.d = obj.classfun(obj.nd, obj.xDisc, obj.d);
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
            
            figure
            hold on
            plot3(obj.x, obj.y, obj.z)
            plot3(obj.x', obj.y', obj.z')
            axis equal
            hold off
            title(sprintf(...
                "NC = [%.2f, %.2f], ND = [%.2f, %.2f]", obj.nc, obj.nd));
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
            a(2) = struct('name', "c", 'min', [1 1], 'max', [10 10]);
            a(3) = struct('name', "d", 'min', 1, 'max', 10);
            a(4) = struct('name', "nc", 'min', [0 0], 'max', [1 1]);
            a(5) = struct('name', "nd", 'min', 0, 'max', 1);
            a(6) = struct('name', "Acu", 'min', [0 0], 'max', [1 1]);
            a(7) = struct('name', "Acl", 'min', [-1 -1], 'max', [0 0]);
            a(8) = struct('name', "Ad", 'min', [0.1 0.1 0.1 0.1], 'max', [1 1 1 1]);            
            
            init_obj = Fuse3d();
            [init_obj, a] = geometry.define(init_obj, a);
        end
        function obj = test()
            
            L = 15;
            c = [10, 10];
            d = 10;
            nc = [0.949, 0.949];
            nd = [0.542, 0.105];
            Acu = flip([0.462, 0.136]);
            Acl = -[0.044, 0.371];
            Ad = [0.258, 0.159, 0.170, 0.199];
            obj = Fuse3d(L, nc, nd, c, d, Acu, Acl, Ad);
            obj.plot()
        end
    end
end


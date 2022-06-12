classdef Fuse3d < Body
    
    properties
        
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
        function self = Fuse3d(length, nc, nd, w, h, Acu, Acl, Ad)
            
            if nargin > 0
                
                self.length = length;
                % Halved as defining half width and upper/lower 
                self.w = w * length/2;
                self.h = h * length/2;
                self.nc = nc;
                self.nd = nd;
                self.Acu = Acu;
                self.Acl = Acl;
                self.Ad = Ad;
            end
            
        end
        function self = generate(self)
            
            S = self.get_shapefuns();
            C = self.get_classfuns();
            
            y_nd = -S.d .* C.d .* (1 - (2 * self.yDisc));
            zu_nd = self.h(1) * C.cu .* S.cu .* (C.d .* S.d);
            zl_nd = self.h(2) * C.cl .* S.cl .* (C.d .* S.d);
            
            % Max or max - min as min will always be 0
            y = y_nd./max(y_nd(:)) * self.w * self.length * 0.5;
            zu = zu_nd./max(zu_nd(:)) * self.h(1) * self.length * 0.5;
            zl = zl_nd./max(-zl_nd(:)) * self.h(end) * self.length * 0.5;
            
            [y, z] = self.combine_upper_lower(zu, zl, y);
            
            x = self.xDisc * self.length;
            x = repmat(x, 1, size(y, 2));
            
            [self.x, self.y, self.z] = self.close_body(x, y, z);
            
            self.upper = zu;
            self.lower = zl;
            self = self.xyz2points();
            self = self.get_data();
        end
        function self = set.Acu(self, in)
            %% TODO: Yay or nae (Change to inc_array fun anyway)
            % Can only be slightly below previous
            for i = 2:numel(in)
                
                if in(i) < in(i-1)*0.9
                    
                    in(i) = in(i-1)*0.9;
                end
            end
            
            self.Acu = in; 
        end
        function self = set.Acl(self, in)
            
            % Can only be slightly above previous
            for i = 2:numel(in)
                
                if in(i) > in(i-1) + 0.1
                    
                    in(i) = in(i-1) + 0.1;
                end
            end
            
            self.Acl = in; 
        end
        function out = interp_shape(self, x, y)
            
            funs = self.shapefuns;
            
            if nargin < 1 || isempty(x)
                
                out(:,1) = funs.Sd;
            else
                out(:,1) = interp1(self.xDisc, funs.Sd, x);
            end
            
            if nargin < 2 || isempty(y)
                
                out(:,2) = funs.Scu;
                out(:,3) = funs.Scl;
            else
                out(:,2) = interp1(self.yDisc, funs.Scu, y);
                out(:,3) = interp1(self.yDisc, funs.Scl, y);
            end
        end
        function out = interp_class(self, x, y)
            
            funs = self.classfuns;
            
            if nargin < 1 || isempty(x)
                
                out(:,1) = funs.Cd;
            else
                out(:,1) = interp1(self.xDisc, funs.Cd, x);
            end
            
            if nargin < 2 || isempty(y)
                
                out(:,2) = funs.Scu;
                out(:,3) = funs.Scl;
            else
                out(:,2) = interp1(self.yDisc, funs.Ccu, y);
                out(:,3) = interp1(self.yDisc, funs.Ccl, y);
            end
        end
        function a = check(self)
            
            upper_diff = self.upper(:,2:end) - self.upper(:,1);
            lower_diff = -(self.lower(:,2:end) - self.lower(:,1));
            
            a = [max(upper_diff(:)), max(lower_diff(:)), ...
                min(self.nose_rad), self.height/self.width];
        end
        function a = get_shapefuns(self)
            
            Au = self.Acu;
            Al = self.Acl;
            
            a.cu = self.shapefun([Au, flip(Au)], self.yDisc);
            a.cl = self.shapefun([Al, flip(Al)], self.yDisc);
            a.d = self.shapefun(self.Ad, self.xDisc);
        end
        function a = get_classfuns(self)
            
            a.cu = self.classfun(self.nc(1), self.yDisc);
            a.cl = self.classfun(self.nc(2), self.yDisc);
            a.d = self.classfun(self.nd, self.xDisc);
        end
        function plotbody(self)
            
            Geometry.plotter(self)
            figure(gcf)
            hold on
            title(sprintf(...
                "NC = [%.2f, %.2f], ND = [%.2f, %.2f]", self.nc, self.nd));
            hold off
        end
        function plotcurves(self)
            
            S = self.get_shapefuns;
            C = self.get_classfuns;
            
            y = -S.d .* C.d;
            zu = C.cu .* S.cu;
            zl = C.cl .* S.cl;
            
%             Geometry.plotter(obj)
            yspare = zeros(size(y));
            zspare = zeros(size(zu));
            figure(gcf)
            hold on
            plot3(self.xDisc, -y, yspare)
            plot3(zspare, self.yDisc-0.5, zu)
            plot3(zspare, self.yDisc-0.5, zl)
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

            a(1) = struct('name', "length", 'min', 10, 'max', 40);
            % Stretch factors (% of length)
            a(2) = struct('name', "w", 'min', 0.1, 'max', 0.8);
            a(3) = struct('name', "h", 'min', [0.1 0.1], 'max', [0.8 0.8]);
            % Spanwise upper/lower curvature: 0 = tall/box, 1 = cylindrical
            a(4) = struct('name', "nc", 'min', [0.1 1e-3], 'max', [1 1]);
            % Streamwise nose/tail curvature: 0 = flat/box, 1 = knife
            a(5) = struct('name', "nd", 'min', [0.25 1e-3], 'max', [0.95 0.8]);
            %% TODO: linear Ac1 > Ac2 rather than fuse height constraint
            a(6) = struct('name', "Acu", 'min', [0.05 0.05], 'max', [1 1]);
            a(7) = struct('name', "Acl", 'min', -[1 1], 'max', -[0.02 0.02]);
            % Have to be > 0 otherwise all(z) == 0
            a(8) = struct('name', "Ad", 'min', [0.3 0.3 0.3 0.3], 'max', [1 1 1 1]);            
            
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
            obj.generate;
            obj.plot()
        end
    end
end


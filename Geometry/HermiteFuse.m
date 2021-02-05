classdef HermiteFuse < Body
    
    properties (SetObservable, AbortSet)
        
        length
        xcp
        ycp
        zcp
        phi
        r
        nradius
        nlength
        noffset = [0 0]
        nose = Nose()
    end
    
    properties
        
        ypoly
        zpoly
    end
    
    properties (Constant)
        
        xPanels = 100
        yPanels = 100
    end
    
    methods
        function obj = HermiteFuse(length, xcp, ycp, zcp, noffset, r, phi)
            
            if nargin > 1
                
                obj.length = length;
                obj.xcp = xcp;
                obj.ycp = ycp;
                obj.zcp = zcp;
                
                if nargin >= 5 && ~isempty(noffset)
                    
                    obj.noffset = noffset;
                end
                
                if nargin >= 6 && ~isempty(r)
                    
                    obj.r = r;
                    obj.phi = phi;
                end
            end
            
            addlistener(obj,'length','PostSet',@obj.initialise);
            addlistener(obj,'xcp','PostSet',@obj.initialise);
            addlistener(obj,'ycp','PostSet',@obj.initialise);
            addlistener(obj,'zcp','PostSet',@obj.initialise);
            addlistener(obj,'noffset','PostSet',@obj.initialise);
            addlistener(obj,'r','PostSet',@obj.initialise);
            addlistener(obj,'phi','PostSet',@obj.initialise);
        end
        function obj = dogenerate(obj)
            
            obj.nose.radius = obj.nradius;
            obj.nose.length = obj.nlength;
            obj.nose.offset = obj.noffset;
            
            [yrads, zrads] = obj.radial_hermite();
            
            dim = size(yrads, 2);
            
            nos = obj.nose;
            nos.yPanels = dim-1;
            
            % Ensure body not interfering with nose
            max_xnose = max(nos.x(end,:)) + 1e-3;
            obj.xcp = max(obj.xcp, max_xnose);
            
            x = [nos.x; repmat(obj.xcp, 1, dim)];
            ymerge = [nos.y; yrads];
            zmerge = [nos.z; zrads];
            
            % Interpolating from end of nose to end of body
            xint = linspace(max_xnose, obj.length, ...
                obj.xPanels)';
            
            for j = dim:-1:1
                
                yint(:,j) = interp1(x(:,j), ymerge(:,j), xint, 'pchip');
                zint(:,j) = interp1(x(:,j), zmerge(:,j), xint, 'pchip');
            end
            
            xint(end+1,:) = xint(end,:) + 1e-9;
            yint(end+1,:) = 0;
            zint(end+1,:) = 0;
            
            obj.x = [nos.x; repmat(xint, 1, dim)];
            obj.y = [nos.y; yint];
            obj.z = [nos.z; zint];
            
            obj.ypoly = yrads;
            obj.zpoly = zrads;
        end
        function [yrads, zrads] = radial_hermite(obj)
            
            row = numel(obj.xcp);
            
            if ~isempty(obj.ycp)
                
                [yrad, zrad] = ...
                    cart2pol(obj.ycp, obj.zcp);
            else
                yrad = repmat(obj.phi, row, 1);
                zrad = obj.r;
            end
            
            int = linspace(0, 1, obj.yPanels);
            
            for i = row:-1:1
                
                while true
                    
                    y = yrad(i,:)';
                    z = zrad(i,:)';
                    
                    mk = finite_diff(y, z);
                    mk(~isfinite(mk)) = 0;
                    mk([1 end],:) = 0;
                    
                    [zp, yp] = hermite(y, z, mk);
                    
                    dim = numel(yp);
                    dim = (0:dim-1)/(dim-1);
                    
                    yp = interp1(dim, yp, int);
                    zp = interp1(dim, zp, int);
                    
                    [yp, zp] = pol2cart(yp, zp);
                  
                    if any(yp < 0)
                        
                        [minr, id] = min(zrad(i, 2:end-1));
                        zrad(i, id+1) = minr + 0.1 * minr;
                    else
                        break
                    end
                end
                
                yrads(i,:) = yp;
                zrads(i,:) = zp;
                
                % test(:,1) = min(y):0.1:max(y);
                % test(:,2) = interp1(y, z, test(:,1), 'pchip');
                % [test(:,1), test(:,2)] = pol2cart(test(:,1), test(:,2));
            end
        end
        function set.xcp(obj, val)
            
            % Normalise val by max
            % Allows redefinition & non-dimensional corrections below
            val = val/max(val);
            
            % Streamwise condition, sets must be 10% away from each other
            for i = 2:numel(val)
                
                comp = val(i-1) + 0.1;
                
                if val(i) < comp
                    
                    val(i) = comp; 
                    
                    % Normalise again if final value was moved
                    if numel(val), val = val/max(val); end
                end
            end
            
            val = obj.length * val(:);
            obj.xcp = val;
        end
        function set.r(obj, val)

            if isvector(val), val = reshape(val, numel(obj.xcp), []); end
            
            [~, col] = size(obj.phi);
            
            % Ensure radially increasing up to max y, decreasing after
            %% TODO: Make generic function
            [~, min_phi_id] = min(abs(obj.phi));
            anchor = val(:,min_phi_id);
            val = min(val, anchor);
            
            for i = 2:col
                
                if i < min_phi_id
                    
                    set = val(:,i) < val(:,i-1);
                    
                elseif i > min_phi_id
                    
                    set = val(:,i) > val(:,i-1);
                else
                    continue
                end
                
                val(set, i) = val(set ,i-1);
            end
            
            %% TODO: Streamwise constraint (> up to mid, < after, or?)
            %% TODO: Allow to be lower than x per metre?
            comp = val(1,:);
            
            for i = 2:size(val, 1)-1
                
                set = val(i,:) < comp;
                val(i, set) = comp(set);
            end
            
            obj.r = obj.length * val;
        end
        function set.ycp(obj, val)
            
            if isvector(val), val = reshape(val, numel(obj.xcp), []); end
            
            % Ensure y increasing radially out to respective max,
            % decreasing after
            [row, col] = size(val);
            [~, maxID] = max(val, [], 2);
            
            for i = 1:row
                for j = 2:col - 1
                    
                    if j < maxID(i)
                        
                        comp = val(i,j-1) + 1e-3;
                        
                        if val(i,j) < comp
                            
                            val(i,j) = comp;
                        end
                        
                    elseif j > maxID(i)
                        
                        comp = val(i,j-1) - 1e-3;
                        
                        if val(i,j) > comp
                            
                            val(i,j) = comp;
                        end
                    end
                end
            end
            
            obj.ycp = val;
        end
        function set.zcp(obj, val)
            
            if isvector(val), val = reshape(val, numel(obj.xcp), []); end
            
            % Ensure z is constantly decreasing
            [~, rad_dim] = size(val);
            
            for i = 2:rad_dim - 1
                
                set = val(:,i) >= val(:,i-1);
                val(set, i) = val(set, i-1) - 1e-3;
            end
            
            obj.zcp = val;
        end
        function set.phi(obj, val)
            
            % Ensure phi is constantly increasing
            for i = 2:numel(val) - 1
                
                comp = val(i-1) - deg2rad(5);
                
                if val(i) > comp, val(i) = comp; end
            end
            
            obj.phi = val;
        end
        function set.noffset(obj, val)
            
            if isempty(obj.zcp)
                
                obj.noffset = -obj.r(1,end) * val;
            else
                obj.noffset = obj.zcp(1, end) * val;
            end
        end
        function doplot(obj)
            
            figure(gcf)
            hold on
            if isempty(obj.ycp)
                
                [y, z] = pol2cart(obj.phi, obj.r);
            else
                y = obj.ycp;
                z = obj.zcp;
            end
            plot3((obj.xcp.*ones(size(obj.ypoly)))', obj.ypoly', obj.zpoly')
            plot3(obj.xcp.*ones(size(y)), y, z, 'kx')
            axis equal
            hold off
        end
    end
    
    methods (Static)
        function [init_obj, a] = define()
            
            ymin = [0 0.05 0.1 0.05 0;
                0 0.05 0.1 0.05 0;
                0 0.05 0.1 0.05 0];
            
            ymax = [0 1 1 1 0;
                0 1 1 1 0;
                0 1 1 1 0];
            
            zmin = [0.1 0.1 -0.5 -1 -1;
                0.1 0.1 -0.5 -1 -1;
                0.1 0.1 -0.5 -1 -1];
            
            zmax = [1 1 0.5 -0.1 -0.1;
                1 1 0.5 -0.1 -0.1;
                1 1 0.5 -0.1 -0.1];
            
            rmin = repmat(0.1, 3, 5);
            rmax = repmat([0.6 0.7 1 0.7 0.6], 3, 1);
            
            pmin = deg2rad([90 5 -40 -85 -90]);
            pmax = deg2rad([90 85 40 -5 -90]);
            
            a(1) = struct('name', "length", 'min', 1, 'max', 20);
            a(2) = struct('name', "xcp", 'min', [0.1 0.4 1], 'max', [0.6, 0.8, 1]);
            a(3) = struct('name', "phi", 'min', pmin(:)', 'max', pmax(:)');
            a(4) = struct('name', "r", 'min', rmin(:)', 'max', rmax(:)');
            a(5) = struct('name', "nradius", 'min', 0, 'max', 0.5);
            a(6) = struct('name', "nlength", 'min', 0.1, 'max', 1);
            a(7) = struct('name', "noffset", 'min', [0 0], 'max', [0 0.9]);
            
            init_obj = HermiteFuse();
            [init_obj, a] = Geometry.define(init_obj, a);
        end
        function obj = test()
            
            leng = 12;
            xcp = [0.3, 0.6, 1]';
            ycp = repmat([0 0.4 0.5 0.41 0], 3, 1);
            zcp = [0.3 0.21 -0.1 -0.21 -0.3;
                0.4 0.21 -0.1 -0.21 -0.4;
                0.3 0.21 -0.1 -0.21 -0.3];
            
            obj = HermiteFuse(leng, xcp, ycp, zcp);
            obj.nradius = 0.25;
            obj.nlength = 0.5;
            obj.noffset(2) = 0.5;
            obj.dogenerate();
            obj.plot;
        end
        function obj = test_polar()
            
            leng = 12;
            xcp = [0.3, 0.6, 1]';
            r = repmat([0.5 0.5 0.8 0.7 0.5], 3, 1);
            phi = deg2rad([90 45 0 -20 -90]);
            
            obj = HermiteFuse(leng, xcp, [], [], [], r, phi);
            obj.nradius = 0.25;
            obj.nlength = 0.5;
            obj.noffset(2) = 0.5;
            obj.dogenerate();
            obj.plot;
            obj.doplot;
        end
        function obj = x34()
            
            leng = 16.431;
            xcp = [4.4238, 16.431]';
            
            ycp = [0 .8793 .91 .86 0];
            zcp = [.8277 -.0187 -.778 -.828 -.828];
            
            j = 1;
            for i = [2, 4]
                
                ycp_int(j,:) = interp1([0 1], ycp(i:i+1), [0.05, 0.5, 0.95]);
                zcp_int(j,:) = interp1([0 1], zcp(i:i+1), [0.05, 0.5, 0.95]);
                j = j + 1;
            end
            
            ycp = [ycp, ycp_int(:)'];
            zcp = [zcp, zcp_int(:)'];
            
            [theta,~] = cart2pol(ycp, zcp);
            [~,id] = sort(theta, 'descend');
            
            ycp = ycp(id);
            zcp = zcp(id);
            
            ycp = repmat(ycp, 2, 1);
            zcp = repmat(zcp, 2, 1);
            
            obj = HermiteFuse(leng, xcp, ycp, zcp);
            obj.nradius = 0.01;
            obj.nlength = 0.1115;
            obj.noffset(2) = .9;
            obj.dogenerate();
            obj.plot;
            obj.doplot;
        end
    end
end
classdef Wingsection < Geometry
    
    properties %(SetObservable)
        
        chord
        span
        trail_sweep
        sections = Aerofoil.get_aerofoil("NACA0006")
        dihedral
        offset = [0 0]
    end
    
    properties
        
        xDiscMethod (1,1) string {mustBeMember(xDiscMethod,{'Linear', 'Cosine', 'Halfcosine'})} = 'Linear'
        xDisc = 100
        % Set to positive number to discretise wing by target dimension
        yDisc = 0
        nParts
        nPanels
    end
    
    properties (Dependent)
        
        area
        taper
        cbar
        lead_edge
        lead_sweep
        upper
        lower
        AR
        Sref
    end
    
    properties (Constant)
        
        name = "wing"
        conical = false
    end
    
    methods
        
        function self = Wingsection(chord, span, sweep, di, sections)
            
            if nargin > 0
                
                self.chord = chord;
                self.trail_sweep = sweep;
                self.dihedral = di;
                self.sections = sections;
                self.span = span;
            end
            
        end
        function self = dogenerate(self)
            
            aerofoils = self.sections;
            nFoil = numel(aerofoils);
            xnorm = self.disc_method();
            
            for i = nFoil:-1:1
                
                aerofoils(i) = aerofoils(i).redist(xnorm);
            end
            
            LE = self.lead_edge;
            
            chords = self.chord;
            nChords = numel(chords);
            
            % Distributing dihedral from partition based to chord based
            di = self.dihedral;
            for i = nChords:-1:1
                if i == 1
                    temp(i) = di(i);
                    
                elseif i == nChords
                    
                    temp(i) = di(end);
                else
                    temp(i) = mean(di(i-1:i));
                end
            end
            di = temp;
            
            if nFoil < nChords
                
                aerofoils(nFoil+1:nChords) = aerofoils(end);
            end
            
            self.nParts = nChords - 1;
            
            xu = [aerofoils.xu] .* chords + LE(:,:,1);
            zu = [aerofoils.zu] .* chords;
            
            xl = [aerofoils.xl] .* chords + LE(:,:,1);
            zl = [aerofoils.zl] .* chords;
            
            yu_rot = zu .* -sin(di) + LE(:,:,2);
            yl_rot = zl .* -sin(di) + LE(:,:,2);
            
            zu_rot = zu .* cos(di) + LE(:,:,3);
            zl_rot = zl .* cos(di) + LE(:,:,3);
            
            self.x = self.wrap(xu, xl) + self.offset(1);
            self.y = self.wrap(yu_rot, yl_rot);
            self.z = self.wrap(zu_rot, zl_rot) + self.offset(2);
            
            if self.yDisc
                
                self = self.span_disc;
            else
                self.nPanels = ones(1, self.nParts);
            end
        end
        function a = check(self, id)
            %% TODO: This (id) or docheck as in Aerofoil
            maxTE = max(self.upper(end,:,3) - self.lower(end,:,3));
            
            % Trying sum to allow for a bit of flexibility
            % jointDiff = min(self.upper(:,1,2) - self.lower(:,1,2));
            jointDiff = mean(self.upper(:,1,2) - self.lower(:,1,2));
            
            dt_ds = diff(max([self.sections.thickness]), [], 2);
            
            a = [min(self.lead_sweep), jointDiff, max(dt_ds), maxTE];
            
            if nargin >= 2 && ~isempty(id)
                
                a = a(id);
            end
        end
        function stag = stagnation(self, U)
            
            unorm = self.quad_data.unit_norm;
            mid = ceil(size(unorm, 2)/2);
            
            del = asin(dotmat(-U, unorm(1,:,:)));
            
            up = del(:,1:mid-1);
            lo = del(:,end:-1:mid+1);
            avg = (up + lo)/2;
            
            [~,id] = max([up; lo; avg]);
            
            id = mode(id);
        end
        function self = set.chord(self, in)
            
            maxTaper = 0.8;
            % Ensuring taper ratio <= 0.8
            for i = 2:numel(in)
                
                if in(i) > in(i-1) * maxTaper
                    
                    in(i) = in(i-1) * maxTaper;
                end
            end
            
            self.chord = in;
        end
        function self = set.trail_sweep(self, in)
            
            deltaMax = deg2rad(20);
            % Constraining sweep on either side of max difference
            for i = 2:numel(in)
                
                in(i) = max(min(in(i), in(i-1) + deltaMax), in(i-1) - deltaMax);    
            end
            
            self.trail_sweep = in;
        end
        function self = set.dihedral(self, in)
            % Dihedral can only increase outboard by max of below
            deltaMax = deg2rad(20);
            
            for i = 2:numel(in)
                
                in(i) = min(max(in(i-1:i)), in(i-1) + deltaMax);
            end
            
            self.dihedral = in;
        end
        function self = set.sections(self, in)
            
            self.sections = self.clean_sections(in);
        end
        function a = get.lead_edge(self)
            
            xTE = [0 cumsum(self.span .* tan(self.trail_sweep))] + self.chord(1);
            yTE = [0 cumsum(self.span .* cos(self.dihedral))];
            zTE = [0 cumsum(self.span .* sin(self.dihedral))];
            
            a(:,:,1) = xTE - self.chord;
            a(:,:,2) = yTE;
            a(:,:,3) = zTE;
        end
        function a = get.lead_sweep(self)
            %% Leading edge sweep calculation
            
            x_le = self.lead_edge(:,:,1);
            diff_le = diff(x_le);
            hyp = (diff_le.^2 + self.span.^2).^0.5;
            
            neg = diff_le < 0;
            
            a = acos(self.span./hyp);
            
            a(neg) = -a(neg);
            a(diff_le == 0) = 0;
        end
        function self = span_disc(self)
            
            target = self.yDisc;
            nPan = 0;
            lo = fun(self.lower);
            up = fun(self.upper);
            
            self.x = self.wrap(up(:,:,1), lo(:,:,1));
            self.y = self.wrap(up(:,:,2), lo(:,:,2));
            self.z = self.wrap(up(:,:,3), lo(:,:,3));
            self.nPanels = nPan;
            
            function b = fun(a)
                
                yz = magmat(a(:,:,[2 3]));
                
                [xDim, yDim] = size(yz);
                
                yzBar = mean(yz, 1);
                dist = diff(yzBar);
                
                % Number of panels, defined once for upper/lower
                % consistency
                if nPan == 0
                    
                    nPan = ceil(dist./target);
                end
                
                % + 1 to give number of points
                b = zeros(xDim, sum(nPan) + 1, 3);
                
                cols = [0 cumsum(nPan)] + 1;
                
                for i = 1:yDim - 1
                    
                    c = cols(i):cols(i+1);
                    b(:,c,:) = a(:,i,:) + (a(:,i+1,:) - a(:,i,:)).*((0:nPan(i))/nPan(i));
                end
            end
        end
        function x_nd = disc_method(self, np, method)
            
            if nargin < 2 || isempty(np), np = self.xDisc; end
            if nargin < 3 || isempty(method), method = self.xDiscMethod; end
            
            x = (0:np)';
            
            if strcmpi(method, "linear")
                
                x_nd = x/max(x);
                
            elseif strcmpi(method, "cosine")
                
                x_nd = 0.5*(1-cos((x*pi)/max(x)));
                
            elseif strcmpi(method, "halfcosine")
                
                x_nd = 1-cos((x*(pi/2))/max(x));
            end
        end
        function self = extendWithin(self)
            % Method to extend wing through to root
            self.x(:,1) = vec_interp(...
                self.y(:,1), ...
                self.y(:,2), ...
                self.x(:,1), ...
                self.x(:,2), ...
                0);
            
            self.z(:,1) = vec_interp(...
                self.y(:,1), ...
                self.y(:,2), ...
                self.z(:,1), ...
                self.z(:,2), ...
                0);
            
            self.x(:,end) = vec_interp(...
                self.y(:,end-1), ...
                self.y(:,end), ...
                self.x(:,end-1), ...
                self.x(:,end), ...
                0);
            
            self.z(:,end) = vec_interp(...
                self.y(:,end-1), ...
                self.y(:,end), ...
                self.z(:,end-1), ...
                self.z(:,end), ...
                0);
            
            self.y(:,1) = 0;
            self.y(:,end) = 0;
            
            self.chord(1) = diff(self.x([1, end], 1));
            self.span(1) = self.y(1,2);
            
%             Geometry.plotter(self)
%             figure(gcf)
%             hold on
%             plot3(x1, y, z1)
%             plot3(xend, y, zend)
        end
        function a = get.upper(self)
            
            p = self.points;
            
            [~,dim,~] = size(p);
            
            a = p(:,1:dim/2,:);
        end
        function a = get.lower(self)
            
            p = self.points;
            
            [~,dim,~] = size(p);
            
            a = fliplr(p(:,(dim/2)+1:end,:));
        end
        function a = get.area(self)
            
            for i = numel(self.span):-1:1
                
                a(i) = 0.5 * sum(self.chord(i:i+1)) * ... 
                    self.span(i) * cos(self.dihedral(i));
            end
        end
        function a = get.taper(self)
            
            a = self.chord(2:end)./self.chord(1:end-1);
        end
        function a = get.cbar(self)
            
            t = self.taper;
            a = (2/3) * self.chord(1:end-1) .* ((1 + t + (t.^2))./(1 + t));
        end
        function a = get.AR(self)
            %% TODO: Remove area = 0 partitions
            a = self.span.^2./self.area;
            a(isnan(a)) = 0;
        end
        function sweep = get_sweep(self, loc)
            
            sweep = atan(tan(self.trail_sweep) - 4./self.AR * ...
                (loc - 1) .* (1 - self.taper)./(1 + self.taper));
            sweep(isnan(sweep)) = 0;
        end
        function section = clean_sections(self, section)
            
            if isobject(section)
                
                dim = self.nParts + 1;
                
                if length(section) ~= dim
                    
                    section(1:dim) = section;
                end
                
            elseif iscell(section)
                
                for i = numel(section):-1:1
                    
                    temp(i) = Aerofoil(section{i});
                end
                
                section = temp;
            end
        end
    end
    
    methods (Static)
        function [init_obj, a] = define()
            
            a(1) = struct('name', "chord", 'min', [0.25 0.2 0.1], 'max', [0.9 0.75 0.7]);
            a(2) = struct('name', "span", 'min', [0.1 0.1], 'max', [0.9 0.7]);
            % a(2) = struct('name', "span", 'min', [0.5 0.5], 'max', [5 5]);
            a(3) = struct('name', "trail_sweep", 'min', -[pi/4 pi/4], 'max', [pi/4 pi/4]);
            a(4) = struct('name', "dihedral", 'min', [0 0], 'max', [pi/4 pi/4]);
            % Scaled by body length
            a(5) = struct('name', "offset", 'min', 0.1, 'max', 0.75);
            
            init_obj = Wingsection();
            [init_obj, a] = Geometry.define(init_obj, a);
        end
        function wing = test()
            
            chord = [5 4 3 1];
            span = [0.5 2 0.5];
            sweep = [0.2 0.4 0];
            di = [0.05 0.1 0];
            
            sections = BezierFoil.test();
            sections(2:length(chord)) = sections;
            
            wing = Wingsection(chord, span, sweep, di, sections);
        end
        function wrapped = wrap(up, lo)
            
            wrapped = [up, fliplr(lo)];
        end
    end
end

classdef Wingsection < Geometry
    
    properties (SetObservable)
        
        chord
        span
        trail_sweep
        sections = Aerofoil.getaerofoil("NACA0006")
        dihedral
        offset = [0 0]
    end
    
    properties
        
        xDisc (1,1) string {mustBeMember(xDisc,{'Linear', 'Cosine', 'Halfcosine'})} = 'Linear'
        % Set to positive number to discretise wing by target dim yDisc
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
        
        function obj = Wingsection(chord, span, sweep, di, sections)
            
            if nargin > 0
                
                obj.chord = chord;
                obj.trail_sweep = sweep;
                obj.dihedral = di;
                obj.sections = sections;
                obj.span = span;
            end
            
            addlistener(obj,'chord','PostSet',@obj.initialise);
            addlistener(obj,'span','PostSet',@obj.initialise);
            addlistener(obj,'trail_sweep','PostSet',@obj.initialise);
            addlistener(obj,'dihedral','PostSet',@obj.initialise);
            addlistener(obj,'sections','PostSet',@obj.initialise);
            addlistener(obj,'offset','PostSet',@obj.initialise);
        end
        function obj = dogenerate(obj)
            
            chords = obj.chord;
            % No rotation applied to wing root
            di = [0 obj.dihedral];
            LE = obj.lead_edge;
            aerofoils = obj.sections;
            
            one = max(length(obj.span), length(obj.trail_sweep));
            two = max(length(obj.chord), length(obj.sections)) - 1;
            
            obj.nParts = max(one, two);
            
            xu = [aerofoils.xu] .* chords + LE(:,:,1);
            zu = [aerofoils.zu] .* chords;
            
            xl = [aerofoils.xl] .* chords + LE(:,:,1);
            zl = [aerofoils.zl] .* chords;
            
            yu_rot = zu .* -sin(di) + LE(:,:,2);
            yl_rot = zl .* -sin(di) + LE(:,:,2);
            
            zu_rot = zu .* cos(di) + LE(:,:,3);
            zl_rot = zl .* cos(di) + LE(:,:,3);
            
            obj.x = obj.wrap(xu, xl) + obj.offset(1);
            obj.y = obj.wrap(yu_rot, yl_rot);
            obj.z = obj.wrap(zu_rot, zl_rot) + obj.offset(2);
            
            if obj.yDisc
                
                obj.span_disc;
            else
                obj.nPanels = ones(1, obj.nParts);
            end
        end
        function a = check(obj, id)
            %% TODO: This (id) or docheck as in Aerofoil
            maxSweep = rad2deg(max(abs(diff(obj.trail_sweep))));
            maxTE = max(obj.upper(end,:,3) - obj.lower(end,:,3));
            
            a = [maxSweep, maxTE];
            
            if nargin >= 2 && ~isempty(id)
                
                a = a(id);
            end
        end
        function stag = stagnation(obj, U)
            
            unorm = obj.quad_data.unit_norm;
            mid = ceil(size(unorm, 2)/2);
            
            del = asin(dotmat(-U, unorm(1,:,:)));
            
            up = del(:,1:mid-1);
            lo = del(:,end:-1:mid+1);
            avg = (up + lo)/2;
            
            [~,id] = max([up; lo; avg]);
            
            id = mode(id);
        end
        function set.chord(obj, in)
            
            % Ensuring taper ratio <= 1
            for i = 2:numel(in)
                
                if in(i) > in(i-1)
                    
                    in(i) = in(i-1);
                end
            end
            
            obj.chord = in;
        end
        function set.sections(obj, in)
            
            obj.sections = obj.clean_sections(in);
        end
        function a = get.lead_edge(obj)
            
            xTE = [0 cumsum(obj.span .* tan(obj.trail_sweep))] + obj.chord(1);
            yTE = [0 cumsum(obj.span .* cos(obj.dihedral))];
            zTE = [0 cumsum(obj.span .* sin(obj.dihedral))];
            
            a(:,:,1) = xTE - obj.chord;
            a(:,:,2) = yTE;
            a(:,:,3) = zTE;
        end
        function a = get.lead_sweep(obj)
            %% Leading edge sweep calculation
            
            x_le = obj.lead_edge(:,:,1);
            diff_le = diff(x_le);
            hyp = (diff_le.^2 + obj.span.^2).^0.5;
            
            neg = diff_le < 0;
            
            a = acos(obj.span./hyp);
            
            a(neg) = -a(neg);
            a(diff_le == 0) = 0;
        end
        function obj = span_disc(obj)
            
            target = obj.yDisc;
            nPan = 0;
            lo = fun(obj.lower);
            up = fun(obj.upper);
            
            obj.x = obj.wrap(up(:,:,1), lo(:,:,1));
            obj.y = obj.wrap(up(:,:,2), lo(:,:,2));
            obj.z = obj.wrap(up(:,:,3), lo(:,:,3));
            obj.nPanels = nPan;
            
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
        function x_nd = disc_method(obj, np, method)
            
            if nargin < 2 || isempty(np), np = 100; end
            if nargin < 3 || isempty(method), method = obj.xDisc; end
            
            x = (0:np)';
            
            if strcmpi(method, "linear")
                
                x_nd = x/max(x);
                
            elseif strcmpi(method, "cosine")
                
                x_nd = 0.5*(1-cos((x*pi)/max(x)));
                
            elseif strcmpi(method, "halfcosine")
                
                x_nd = 1-cos((x*(pi/2))/max(x));
            end
        end
        function a = get.upper(obj)
            
            p = obj.points;
            
            [~,dim,~] = size(p);
            
            a = p(:,1:dim/2,:);
        end
        function a = get.lower(obj)
            
            p = obj.points;
            
            [~,dim,~] = size(p);
            
            a = fliplr(p(:,(dim/2)+1:end,:));
        end
        function a = get.area(obj)
            
            for i = numel(obj.span):-1:1
                
                a(i) = 0.5 * sum(obj.chord(i:i+1)) * ... 
                    obj.span(i) * cos(obj.dihedral(i));
            end
        end
        function a = get.taper(obj)
            
            a = obj.chord(2:end)./obj.chord(1:end-1);
        end
        function a = get.cbar(obj)
            
            t = obj.taper;
            a = (2/3) * obj.chord(1:end-1) .* ((1 + t + (t.^2))./(1 + t));
        end
        function a = get.AR(obj)
            %% TODO: Remove area = 0 partitions
            a = obj.span.^2./obj.area;
            a(isnan(a)) = 0;
        end
        function sweep = get_sweep(obj, loc)
            
            sweep = atan(tan(obj.trail_sweep) - 4./obj.AR * ...
                (loc - 1) .* (1 - obj.taper)./(1 + obj.taper));
            sweep(isnan(sweep)) = 0;
        end
        function section = clean_sections(obj, in)
            
            if isobject(in)
                
                dim = obj.nParts + 1;
                
                if length(in) ~= dim
                    
                    in(1:dim) = in;
                end
                
            elseif iscell(in)
                
                for i = numel(in):-1:1
                    
                    temp(i) = Aerofoil(in{i});
                end
                
                in = temp;
            end
            
            xnorm = obj.disc_method();
            for i = numel(in):-1:1
                
                section(i) = in(i).redist(xnorm);
            end
        end
    end
    
    methods (Static)
%         function plot()
%             
%         end
        function [init_obj, a] = define()
            
            a(1) = struct('name', "chord", 'min', [0.4 0.25 0.1], 'max', [1.5 1 0.75]);
            a(2) = struct('name', "span", 'min', [1 0.1], 'max', [2 2]);
            % a(2) = struct('name', "span", 'min', [0.5 0.5], 'max', [5 5]);
            a(3) = struct('name', "trail_sweep", 'min', [0 0], 'max', [1.4 1.4]);
            a(4) = struct('name', "dihedral", 'min', [0 0], 'max', [pi/4 pi/4]);
            a(5) = struct('name', "offset", 'min', [0.1 -0.4], 'max', [0.7 0.4]);
            
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

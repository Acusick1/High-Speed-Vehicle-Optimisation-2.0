classdef Bezier < handle
    
    properties (SetObservable, AbortSet)
        
        control_points
    end
    properties
        
        order
        
        %% TODO: Linear constraints solved by optimiser, not geometries
        con_xcp = true
        % Does nothing. Implemented as optimiser constraint instead
        con_zend = true
    end
    properties (Constant)
        
        s = (0:0.01:1)';
    end
    properties (Dependent)
        
        xcp
        ycp
        curve
    end
    
    methods
        function obj = Bezier(varargin)
            
            if nargin > 0
                
                if nargin == 1
                    
                    obj.control_points = varargin{1};
                else
                    dim = length(varargin);
                    obj = repelem(obj, dim);
                    
                    for i = 1:dim
                        
                        obj(i) = obj.array(varargin{i});
                    end
                end
            end
            
            %% For testing
            % addlistener(obj,'control_points','PostSet',@obj.update_plot);
        end
        function set.control_points(obj, val)
            % Input: Vector of [xb1, zb1, ... xbn, zbn]
            %        Matrix of [xb11, zb11, ... xbn1, zbn1
            %                               ...
            %                   xb1m, zb1m, ... xbnm, xbnm]
            
            if isvector(val)
                
                val = reshape(val, obj.order + 1, []);
            end
            
            %% TODO: Constraint hacks as in properties
            % Every second row is x poisition, constrain if required
            if obj.con_xcp
                
                next_row = val(end, 1:2:end);
                
                for i = size(val, 1)-1:-1:1
                    
                    row = val(i, 1:2:end);
                    con = row < next_row;
                    row(con) = next_row(con);
                    val(i, 1:2:end) = row;
                    
                    next_row = row;
                end
            end
            
            obj.control_points = val;
        end
        function a = get.xcp(obj)
            
            a = obj.control_points(:,1:2:end);
        end
        function a = get.ycp(obj)
            
            a = obj.control_points(:,2:2:end);
        end
        function a = get.order(obj)
            
            a = size(obj.xcp, 1) - 1;
        end
        function a = get.curve(obj)
            
            x = obj.xcp;
            y = obj.ycp;
            
            n = obj.order;
            k = 0:n;
            % Binomial coefficients
            f = factorial(n)./(factorial(n-k) .* factorial(k));
            % Component shape functions
            b = f.*((1 - obj.s).^(n - k)).*(obj.s.^k);
            
            for i = size(x, 2):-1:1
                
                a(:,1,i) = flipud(sum(b .* x(:,i)', 2));
                a(:,2,i) = flipud(sum(b .* y(:,i)', 2));
            end
        end
        function plot(obj)
            
            figure(gcf)
            hold on
            h = plot(obj.xcp, obj.ycp, 'o');
            
            for i = 1:length(h)
                
                plot(obj.curve(:,1,i), obj.curve(:,2,i), 'color', h(i).Color);
            end
            hold off
        end
        function update_plot(obj, varargin)
            
            if ishandle(1)
                
                figure(1)
                hold on
                clf
                obj.plot
            end
        end
    end
    methods (Static)
        function obj = array(in)
            
            obj = Bezier(in);
        end
    end
end
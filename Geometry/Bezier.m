classdef Bezier < handle
    
    properties (SetObservable, AbortSet)
        
        control_points
    end
    properties
        
        order
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
        end
        function set.control_points(obj, val)
            
            if isvector(val)
                
                val = reshape(val, obj.order, []);
            end
            
            obj.control_points = val;
        end
        function a = get.xcp(obj)
            
            a = obj.control_points(:,1:2:end);
        end
        function a = get.ycp(obj)
            
            a = obj.control_points(:,2:2:end);
        end
        function a = get.curve(obj)
            
            x = obj.xcp;
            y = obj.ycp;
            
            %% TODO: CHECK WHAT ORDER SHOULD BE (ncp | - 1)
            n = size(obj.xcp, 1);
            
            k = 0:n - 1;
            f = factorial(n)./(factorial(k).*factorial(n-k));
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
    end
    methods (Static)
        function obj = array(in)
            
            obj = Bezier(in);
        end
    end
end
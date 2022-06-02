classdef Flightstate
    
    properties
        
        alpha
        Minf
        altitude
        beta = 0
        delta
        Mach_q
        delta_q
        Pr
        Uinf
        Tinf
        rinf
        mu
        kt
        a
        Pinf
        q
        order (1,1) string {mustBeMember(order,{'full_fact', 'None'})} = 'full_fact'
        max_beta
        max_delta
    end
    
    properties (Constant)
        
        gamma = 1.4;
        R = 287;
        cp = Flightstate.R * Flightstate.gamma/(Flightstate.gamma-1)
    end
    
    properties (Dependent)

        U
    end
    
    methods
        function obj = Flightstate(varargin)
            %% Flightstate class
            % Inputs: alpha, Mach, altitude
            
            if nargin > 0
                
                if strcmpi(obj.order, "full_fact")
                
                    states = Flightstate.full_fact(varargin);
                else
                    for i = numel(varargin):-1:1
                    
                        states(:,i) = 1:numel(varargin);
                    end
                end
                
                [~, table] = theta_beta_mach_curves();
                fn = fieldnames(Flightstate);
                
                for i = size(states, 1):-1:1
                    for j = 1:size(states, 2)
                        
                        obj(i).(fn{j}) = varargin{j}(states(i,j));
                    end
                    
                    obj(i) = obj(i).get_atmospheric_values();
                    obj(i) = obj(i).max_shock_angles(table);
                end
            end
        end
        
        function obj = get_atmospheric_values(obj)
                
            output = tewari_atmosphere(obj.altitude, 0, 0);
            
            obj.Tinf = output(1); % Freestream temperature
            obj.rinf = output(2); % Freestream density
            obj.Pinf = output(7); % Freestream pressure
            obj.a = output(5); % Speed of sound
            obj.mu = output(8); % Dynamic viscosity
            obj.kt = output(10); % Thermal conductivity
            
            g = obj.gamma;
            M = obj.Minf;
            
            % Freestream stagnation pressure ratio
            Pinf_P0 = (2./((g+1)*(M.^2))).^(g/(g-1)) .* ...
                (((2*g*(M.^2))-(g-1))/(g+1)).^(1/(g-1));
            
            % Matching points for Newtonian + Prandtl-Meyer method
            % CHECK: reasonable results for all Mach numbers?
                
            [obj.delta_q, obj.Mach_q] = matching_point(g, Pinf_P0);
            
            obj.Pr = obj.mu .* obj.cp./obj.kt; % Prandtl number
            obj.Uinf = M .* obj.a;
            obj.q = 0.5 * obj.rinf .* obj.Uinf.^2;
        end
        function obj = max_shock_angles(obj, table)
            
            if nargin < 1
             
                [~, table] = theta_beta_mach_curves();
            end
            
            [~, angles] = halfspace(obj.Minf, table);
            
            obj.max_delta = angles(:,1);
            obj.max_beta = angles(:,2);
        end
        function a = get.U(obj)
            
            a = obj.Uinf * [cos(obj.alpha) * cos(obj.beta), sin(obj.beta), sin(obj.alpha)];
        end
    end
    
    methods (Static)
        function states = full_fact(def)
            
            for i = numel(def):-1:1
                
                nVar(i) = numel(def{i});
            end
            
            nVar = nVar(nVar > 0);
            states = fullfact(nVar);
        end
    end
end
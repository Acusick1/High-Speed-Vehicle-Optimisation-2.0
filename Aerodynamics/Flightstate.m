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
        function self = Flightstate(varargin)
            %% Flightstate class
            % Inputs: alpha, Mach, altitude
            
            if nargin > 0
                
                if strcmpi(self.order, "full_fact")
                
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
                        
                        self(i).(fn{j}) = varargin{j}(states(i,j));
                    end
                    
                    self(i) = self(i).get_atmospheric_values();
                    self(i) = self(i).max_shock_angles(table);
                end
            end
        end
        
        function self = get_atmospheric_values(self)
                
            output = tewari_atmosphere(self.altitude, 0, 0);
             
            self.Tinf = output(1); % Freestream temperature
            self.rinf = output(2); % Freestream density
            self.Pinf = output(7); % Freestream pressure
            self.a = output(5); % Speed of sound
            self.mu = output(8); % Dynamic viscosity
            self.kt = output(10); % Thermal conductivity
            
            g = self.gamma;
            M = self.Minf;
            
            % Freestream stagnation pressure ratio
            Pinf_P0 = (2./((g+1)*(M.^2))).^(g/(g-1)) .* ...
                (((2*g*(M.^2))-(g-1))/(g+1)).^(1/(g-1));
            
            % Matching points for Newtonian + Prandtl-Meyer method
            % CHECK: reasonable results for all Mach numbers?
                
            [self.delta_q, self.Mach_q] = matching_point(g, Pinf_P0);
            
            self.Pr = self.mu .* self.cp./self.kt; % Prandtl number
            self.Uinf = M .* self.a;
            self.q = 0.5 * self.rinf .* self.Uinf.^2;
        end
        function self = max_shock_angles(self, table)
            
            if nargin < 1
             
                [~, table] = theta_beta_mach_curves();
            end
            
            [~, angles] = halfspace(self.Minf, table);
            
            self.max_delta = angles(1);
            self.max_beta = angles(2);
        end
        function a = get.U(self)
            
            a = self.Uinf * [cos(self.alpha) * cos(self.beta), sin(self.beta), sin(self.alpha)];
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
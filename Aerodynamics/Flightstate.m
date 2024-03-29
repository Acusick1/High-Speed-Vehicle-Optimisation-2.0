%FLIGHTSTATE defines a vehicles attitude towards the freestream flow and
%the necessary atmospheric properties for aerodynamic computations
classdef Flightstate
    
    properties
        
        alpha       % Angle of attack (rad)
        Minf        % Freestream Mach number
        altitude    % Altitude (m)
        beta = 0    % Sideslip angle (rad)
        delta       % Control surface deflection angle (rad)
        %% TODO: move to aero? Need both?
        Mach_q      % Matching point Mach number
        delta_q     % Matching point deflection angle
        %%
        Pr          % Prandtl number
        Uinf        % Freestream velocity (m/s)
        Uvec        % Freestream velocity vector [3x1] (m/s) 
        Tinf        % Freestream temperature
        Pinf        % Freestream pressure
        rinf        % Freestream density
        mu          % Dynamic viscosity
        kt          % Thermal conductivity
        a           % Local speed of sound        
        q           % Dynamic pressure
        T0          % Stagnation temperature
        r           % Recovery factor [laminar, turbulent]
        
        % Table of maximum shockwave and deflection angles for range of 
        % Mach numbers
        max_theta_beta_mach
        
        max_beta    % Max weak shockwave angle for freestream Mach number
        max_del     % Max deflection angle for freestream Mach number
    end
    
    properties (Constant)
        
        gamma = 1.4;
        STF = 5.67e-8;  % Stefan-Boltzman constant
        R = 287;
        S = 110         % Sutherland's constant
        eps = 0.8       % Emissivity
        cp = Flightstate.R * Flightstate.gamma/(Flightstate.gamma-1)
    end
    
    methods
        function self = Flightstate(varargin)
            %FLIGHTSTATE constructor
            %   Inputs:
            %   alpha - Angle of attack (rad)
            %   Minf - Freestream Mach number
            %   altitude - Altitude (m)
            %   order - Optional string to define how arrays of inputs are
            %       handled:
            %       "fullfact" (default) will make a Flightstate object for 
            %           every combination of the previously defined inputs.
            %       "none" will make Flightstate objects of the inputs as
            %           given e.g. if each array if length 3, three objects
            %           will be created. Note that inputs with length less
            %           than the length of the maximum input will be have
            %           their final element used, this allows single inputs
            %           to be combined with variable inputs e.g. 3 angle of
            %           attacks at 1 Mach number.
            
            if nargin > 0
                
                % Extract order argument if present
                order_id = cellfun(@isstring, varargin);
                order_arg = varargin(order_id);
                % Ensure arguments are without optional order argument
                args = varargin(~order_id);
                
                % Get length of each argument
                arg_lengths = cellfun(@numel, args);
                
                if any(arg_lengths > 1)
                    % If more than one argument set present, move to
                    % wrapper function
                    self = Flightstate.define_many(...
                        args, arg_lengths, order_arg{:});
                else
                    % Define object with singular inputs
                    self.alpha      = args{1};
                    self.Minf       = args{2};
                    self.altitude   = args{3};
                    
                    if numel(args) >= 4, self.beta = args{4}; end
                    
                    [~, self.max_theta_beta_mach] = theta_beta_mach_curves();
                    self = self.max_shock_angles();
                    self = self.get_atmospheric_values();
                end
            end
        end
        function self = set.Minf(self, value)
            % Set Minf and run dependent functions
            self.Minf = value;
            self = self.max_shock_angles();
            self = self.get_atmospheric_values();
        end
        function self = set.altitude(self, value)
            % Set altitude and run dependent functions
            self.altitude = value;
            self = self.get_atmospheric_values();
        end 
        function self = get_atmospheric_values(self)
            %GET_ATMOSPHERIC_VALUES derives necessary atmospheric
            %properties used in aerodynamic calculations
            
            % Ensure this function does not run unless both mach and 
            % altitude are set
            if isempty(self.Minf) || isempty(self.altitude), return; end
            
            output = tewari_atmosphere(self.altitude, 0, 0);
             
            self.Tinf = output(1);  % Freestream temperature
            self.rinf = output(2);  % Freestream density
            self.Pinf = output(7);  % Freestream pressure
            self.a = output(5);     % Speed of sound
            self.mu = output(8);    % Dynamic viscosity
            self.kt = output(10);   % Thermal conductivity
            
            g = self.gamma;
            M = self.Minf;
            
            % Freestream stagnation pressure ratio
            Pinf_P0 = (2./((g+1)*(M.^2))).^(g/(g-1)) .* ...
                (((2*g*(M.^2))-(g-1))/(g+1)).^(1/(g-1));
            
            % Matching point method for Newtonian + Prandtl-Meyer method
            [self.delta_q, self.Mach_q] = matching_point(g, Pinf_P0);
            
            self.Pr = self.mu .* self.cp./self.kt; % Prandtl number
            % Recovery factor for laminar and turbulent boundary layers
            self.r = [self.Pr^(1/2), self.Pr^(1/3)];
            self.Uinf = M .* self.a;
            self.q = 0.5 * self.rinf .* self.Uinf.^2;
            % Stagnation temperature
            self.T0 = self.Tinf * (1 + (g - 1)/2 * (M^2));
            
        end
        function self = max_shock_angles(self)
            %MAX_SHOCK_ANGLES Defines maximum deflection and shockwave 
            %angles for input Mach number from theta beta Mach curves
            [~, angles] = halfspace(self.Minf, self.max_theta_beta_mach);
            
            self.max_del = angles(1);
            self.max_beta = angles(2);
        end
        function Uvec = get.Uvec(self)
            % Putting this in a get function as dependent on multiple
            % properties, will have to be updated if any are updated
            if ~isempty(self.Uinf)
                
                Uvec = self.Uinf * ...
                    [cos(self.alpha) * cos(self.beta),...
                     sin(self.beta), ...
                     sin(self.alpha)];
            end
        end
    end
    
    methods (Static)
        function obj = define_many(all_inputs, dimensions, order)
            %DEFINE_MANY is a wrapper function to define multiple
            %Flightstate objects from input arrays
            if nargin < 3 || isempty(order), order = "fullfact"; end
            
            if strcmpi(order, "fullfact")
                % Get all combination IDs for each input
                input_id = fullfact(dimensions);
            elseif strcmpi(order, "none")
                % Use inputs as given, extend any inputs as necessary
                max_dim = max(dimensions);
                
                % Defining inputs per flightstate
                % If an input has less elements than input with the most 
                % elements, the final element will be repeated
                for i = numel(all_inputs):-1:1
                    input_id(:,i) = min(1:max_dim, dimensions(i));
                end
            else
                error("No setup for multiple Flightstate order given: %s", order)
            end
            
            % Rows represent set of input IDs for a single flight state
            for i = size(input_id, 1):-1:1
                
                % Indexing each input to create a single flight state set
                for j = size(input_id, 2):-1:1
                    state_inputs{j} = all_inputs{j}(input_id(i,j));
                end
                
                obj(i) = Flightstate(state_inputs{:});
            end
        end
    end
end
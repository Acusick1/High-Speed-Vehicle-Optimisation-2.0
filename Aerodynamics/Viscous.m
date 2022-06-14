classdef Viscous
    
    properties
        
        Minf            % Freestream flow properties
        Pinf
        rinf
        Tinf
        mu
        a               % Speed of sound at altitude
        eps = 0.8       % Emissivity
        Pr              % Prandtl number
        turbulent
        
        part
        del             % Flow inclination
        area
        Rn              % Nose radius
        Aref
        
        Tedge
        Medge
        Vedge
        Pedge
        rho
        Re
        
        Twall
        qwall
        max_Twall
        cf
        
        Cf
        Re_x       % Characteristic length (for Re_x: distance along stream)
        Taw     % Adiabatic wall temperature
    end
    
    properties (Dependent)
        
        T0      % Stagnation temperature
        Uinf
        r       % Recovery factor
        Re_T    % Transition reynolds
    end
    
    properties (Constant)
        
        STF = 5.67e-8;  % Stefan-Boltzman constant
        gamma = 1.4;
        R = 287;
        cp = Flightstate.R * Flightstate.gamma/(Flightstate.gamma-1)
        S = 110         % Sutherland's constant
    end
    
    methods
        function self = Viscous()
            %% TODO: Implement this, set constant values to remove remaining get functions
        end
        function self = test(self)
            
            dim = size(self.Tedge);
            
            if isempty(self.turbulent)
                
                self.turbulent = true(dim);
            end
            % Adiabatic wall condition
            % Tw = self.get_Taw(1,:);
            
            % Cold wall condition
            Tw = zeros(dim);
            if ~isempty(self.Twall)
            
                Tw = Tw + self.Twall;
                self.Twall = Tw;
                max_it = 1;
            else
                max_it = 100;
            end
            set = self.Tedge <= 0;
            
            if ~all(set(:))
                
                self.Tedge(set) = min(self.Tedge(~set));
                self.rho(set) = min(self.rho(~set));
            end
            
            self.Re_x = self.get_Re_x();
            
            % [self.cf, q, self.Re] = obj.eckert(Tw);
            
            for i = 1:max_it
                %% TODO: Comapred to DiGiorgio2019, commenting out gives closer result, cold wall assumption?
                %% TODO: Also producing imaginary numbers
                [self.cf, q, self.Re] = self.eckert(Tw);
                
                prev = Tw;
                % qw = self.simple_heating(Tw, Rn);
                
                % Aerothermodynamics of Transatmospheric Vehicles
                % Aero-Thermodynamics for Conceptual Design
                % Tw = real((q/(self.eps * self.STF)).^0.25);
                
                % DEPLOYED PAYLOAD ANALYSIS FOR A SINGLE STAGE TO ORBIT SPACEPLANE
                % Sizing of Conceptual Hypersonic Long-Range Transport Aircraft using a Multi-Disciplinary Optimisation Strategy
                % Tw = real((self.Tinf^4 + (qw/(self.eps * self.STF))).^0.25);
                
                % Comparison of Engineering Correlations for Predicting Heat Transfer in Zero-pressure-gradient Compressible Boundary Layers with CFD and Experimental Data
                % Tw = ((q - qw)/(self.eps * self.STF) + self.Tinf^4).^0.25;
                Tw = real(((q - 0)./(self.eps * self.STF) + self.Tinf.^4).^0.25);
                
                %% TODO: CLEAN
                Tw(self.part.data.area == 0) = 0;
                %con = Tw > self.get_Taw | ~isfinite(Tw) | imag(Tw) > 0;
                %Tw(con) = self.get_Taw(con);
                
                if i > 1
                    %% TODO: Leading edge instability
                    dTw = Tw(2:end,:) - prev(2:end,:);
                    abs_diff = sum(dTw(:).^2);
                    
                    if all(isfinite(Tw(:))) && (abs_diff < 1e-6 || i >= max_it)
                        
                        % Only reset Twall here incase of prescribed Twall
                        self.Twall = Tw;
                        break
                    end
                end
            end
            
            self.qwall = q;
            self.max_Twall = max(Tw(:));
            self.Cf = self.get_Cf();
        end
        function [cf_star, q, Re_star_x] = eckert(self, Tw)
            
            Te = self.Tedge;
            Me = self.Medge;
            % Pe = self.Pre;
            % Ve = Me * self.a;
            Ve = magmat(self.Vedge);
            turb = self.turbulent;
            
            %%
            % Original Eckert's reference temperature method
            % T_star = Te .* (1 + 0.032*(Me.^2) + 0.58*((Tw./Te) - 1));
            
            % Thermoelastic Formulation of a Hypersonic Vehicle Control Surface for Control-Oriented Simulation
            % Tt = Te .* (1 + (self.gamma - 1) .* ((Me.^2)/2));
            % Tr = self.r(2) * (Tt - Te) + Te;
            % T_star = Te + 0.5*(Tw - Te) + 0.22*(Tr - Te);
            
            % Meador-Smart reference temperature method - Lam & Turb
            T_star1 = Te .* (0.45 + 0.55*(Tw./Te) + 0.16*self.r(1) * ...
                ((self.gamma - 1)/2) .* Me.^2);
            
            T_star2 = Te .* (0.5*(1 + (Tw./Te)) + 0.16*self.r(2) * ...
                ((self.gamma - 1)/2) .* Me.^2);
            
            T_star = zeros(size(turb));
            T_star(~turb) = T_star1(~turb);
            T_star(turb) = T_star2(turb);
            
            % Ideal gas law
            % Pedge Anderson2006 p368
            rho_star = self.Pedge./(self.R * T_star);
            % rho_star = self.Pinf./(self.R * T_star);
            
            % Sutherland's law
            mu_star = self.mu * ((T_star/self.Tinf).^1.5) .* ...
                (self.Tinf + self.S)./(T_star + self.S);
            
            Re_star_x = (rho_star .* Ve .* self.Re_x)./mu_star; 
            Re_star_x(Re_star_x < 1) = 0;
            
            % Accurate up to Re 10^9 according to 
            % Thermoelastic Formulation of a Hypersonic Vehicle Control Surface for Control-Oriented Simulation
            % cf_star = (0.37./(log(Re_star_x).^2.584));
            
            % Equation from above ref
            % cf_star = (0.0592./(Re_star_x.^0.2));
            
            % Hypersonic and high-temperature gas dynamics (Meador-Smart)
            cf_star = zeros(size(turb));
            cf_star(~turb) = 0.664./((Re_star_x(~turb)).^0.5);
            cf_star(turb) = 0.02296./((Re_star_x(turb)).^0.139);
            
            % Mangler factor for conical bodies
            if self.part.conical
                
                % White2006
                ml = 1;
                cf_star(~turb) = (2 + ml)^(ml/(ml+1)) * cf_star(~turb);
                
                mt = 0.25;
                cf_star(turb) = (2 + mt)^(mt/(mt+1)) * cf_star(turb);
            end
            
            % redge = self.rho;
            %% TODO: Another (more reputable?) source also
            % Thermoelastic Formulation of a Hypersonic Vehicle Control Surface for Control-Oriented Simulation
            % Anderson2006 p286 (but stated as qw)
            St_star = 0.5 * cf_star * self.Pr^(-2/3);
            % q = redge .* Ve .* St_star .* (self.get_Taw - Tw) * self.cp;
            q = rho_star .* Ve .* St_star .* (self.get_Taw - Tw) * self.cp;
            
            % Can be nan/inf due to zero area panels etc
            con = ~isfinite(cf_star);
            cf_star(con) = 0;
            q(con) = 0;
        end
        function qw = simple_heating(self, Tw, Rn)
            
            if nargin < 1 || isempty(Tw)
                % Tw = self.Tinf; 
                Tw = 0;
            end
            
            h0 = self.get_enthalpy(self.T0);
            hw = self.get_enthalpy(Tw);
            d = self.del;
            
            [rows, cols] = size(d);
            qw = zeros(rows, cols);
            
            con = self.Re < self.get_Re_T;
            % Ensures below will work as expected (one true per column)
            % Also ensures that in the work case strip will be set as fully
            % turbulent rather than fully laminar
            con(1,:) = true;
            
            % Condition transposed so that first of every column found
            % rather than row, output matched
            [col, row] = find(~con', cols, 'first');
            
            % ids are for first Re_turb, therefore start xt from previous
            % setting first row as true allows this consistently
            row = row - 1;
            
            xt_id = sort((col - 1) * rows + row);
            xt = self.Re_x(xt_id)';
            
            qw_lam = lam();
            qw_turb = turb(self.Re_x - xt);
            
            qw(con) = qw_lam(con);
            qw(~con) = qw_turb(~con);
            qw(~isfinite(qw)) = 0;
            
            if ~isempty(Rn)
                
                q0 = stagnation(Rn);
                % Only invoked if nose radius is positive, and panel faces
                % flow, otherwise flap plate relation holds
                con = Rn > 0 & self.Tedge(1,:) > self.Tinf;
                
                if self.part.conical
                    
                    qw(1,con) = q0(con);
                else
                    %% TODO: Move to Wingsection
                    % Hacky way to get sweep for upper/lower surfaces
                    le_sweep = self.part.get_lead_sweep();
                    le_sweep = le_sweep([1:end end end:-1:1]);
                    
                    % From: Aerothermodynamics of Transatmospheric Vehicles
                    qw(1,con) = (0.5 * q0(con).^2 .* cos(le_sweep(con)).^2 ...
                        + qw(1,con).^2 .* sin(le_sweep(con)).^2).^0.5;
                end
            end
            
            function qw = stagnation(Rn)
                
                M = 3;
                N = 0.5;
                C = 1.83e-8 * Rn.^-0.5 .* (1 - hw(1,:)./h0);
                
                qw = self.rinf^N * self.Uinf^M * C;
            end
            function qw = lam(row)
                
                if nargin == 1 && ~isempty(row)
                    
                    d = d(row,:);
                    X = self.Re_x(row,:);
                else
                    X = self.Re_x;
                end
                
                M = 3.2;
                N = 0.5;
                C = 2.53e-9 * cos(d).^0.5 .* sin(d) .* ...
                    X.^-0.5 .* (1 - hw/h0);
                
                qw = self.rinf^N * self.Uinf^M * C;
            end
            function qw = turb(xt)
                %% TODO: Skin friction hack
                % Using abs(sin(d)) as value is equal but negative when d 
                % negative. This would result in negative skin friction
                % over such panels.
                % However now assuming -d gives same skin friction as +d,
                % At least for this term
                
                if self.Uinf > 3962                    
                
                    M = 3.7;
                    N = 0.8;
                    C = 2.2e-9 * cos(d).^2.08 .* sin(abs(d)).^1.6 .* ...
                        xt.^-0.2 .* (1 - (1.11*hw/h0));
                else
                    % Need finite value otherwise Inf
                    Tw(Tw == 0) = 1e-16;
                    
                    M = 3.37;
                    N = 0.8;
                    C = 3.89e-8 * cos(d).^1.78 .* sin(abs(d)).^1.6 .* ...
                        xt.^-0.2 .* (Tw/566).^-0.25 .* (1 - (1.11*hw/h0));
                end
                
                qw = self.rinf^N * self.Uinf^M * C;
            end
        end
        
        function a = get_Re_T(self)
            
            % Source: Three-Dimensional Laminar Boundary-Layer Transition on a Sharp 88 Cone at Mach 10
            a = exp(6.421 * exp(1.209e-4 * self.Medge.^2.641));
        end
        
        function qw = medium_heating(self, Tw, Rn, laminar)
            
            if nargin < 1 || isempty(Tw), Tw = self.Tinf; end
            if nargin < 3 || isempty(laminar), laminar = true; end
            
            redge = self.Pedge/(287 * self.Tedge);
            
            if ~isempty(Rn)
                
                if self.part.conical
                    
                    qw = stagnation(Rn);
                else
                    % Hacky way to get sweep for upper/lower surfaces
                    le_sweep = self.part.get_lead_sweep();
                    le_sweep = le_sweep([1:end end end:-1:1]);
                    
                    q0 = stagnation(Rn);
                    qfp = lam(1);
                    
                    qw = (0.5 * q0.^2 .* cos(le_sweep).^2 + ...
                        qfp.^2 .* sin(le_sweep).^2).^0.5;
                end
            else
                
            end
            
            function stagnation(Rn)
                
                due_dx = 1./Rn .* (2 * (self.Pedge - self.Pinf)./redge).^0.5;
                
                if self.conical
                    
                    const = 0.763;
                else
                    const = 0.57;
                end
                
                qw = const * self.Pr^-0.6 * (redge * medge).^0.5 ...
                    * due_dx.^0.5 .* (haw - hw);
            end
        end
        function h = get_enthalpy(self, T)
            % Calorifically perfect gas: Anderson2006 p276
            h = self.cp * T;
        end
        function T = get_temperature(self, h)
            % Calorifically perfect gas: Anderson2006 p276
            T = h / self.cp;
        end
        function Cf = get_Cf(self)
            
            part_area = self.part.data.area;
             
            dimensionalise = (self.cf .* part_area)/self.Aref;
            Cf = sum(dimensionalise(:));
        end
        function Re_x = get_Re_x(self)
            %% TODO: Move to Geometry
            if self.part.conical
                
                le = mean(self.part.points(1,:,:));
                centre = self.part.centre - le;
                centre = [zeros(1, size(centre, 2), 3); centre];
            else
                %% TODO: Correct? Make consistent                
                try
                    le = (self.part.points(1,1:end-1,:) + ...
                        self.part.points(1,2:end,:))/2;
                    centre = self.part.centre - le;
                    centre = [zeros(1, size(centre, 2), 3); centre];
                catch
                    le = zeros(1, 1, 2);
                    centre = self.part.data.centre - le;
                    centre = [zeros(1, size(centre, 2), 2); centre];
                end
            end

            Re_x = cumsum(magmat(diff(centre)));
        end
        function a = get.Uinf(self)
            
            a = self.Minf * self.a;
        end
        function a = get.T0(self)
            
            a = self.Tinf * (1 + (self.gamma - 1)/2 * (self.Minf^2));
        end
        function a = get_Taw(self)
            
            con = self.turbulent;
            
            a = zeros(size(con));
            
            a(con) = self.r(2)*(self.T0 - self.Tedge(con)) + self.Tedge(con);
            a(~con) = self.r(1)*(self.T0 - self.Tedge(~con)) + self.Tedge(~con);
        end
        function a = get.r(self)
            
            % [laminar, turbulent]
            a = [self.Pr^(1/2), self.Pr^(1/3)];
        end
    end
    
    methods (Static)
        
        function obj = from_aerodynamics(aero)
            %% Create viscous object from existing aerodynamics object
            % Apply fields from existing Aerodynamics obj to Viscous obj
            
            obj = Viscous();
            aero_fn = fieldnames(aero);
            visc_fn = fieldnames(obj);
            
            for i = 1:numel(visc_fn)
                
                field = visc_fn{i};
                
                if any(strcmp(field, aero_fn))
                    
                    obj.(field) = aero.(field);
                end
            end
            
            %% TODO: Redefine aero/flightstate/viscous to avoid interdep
            flow = aero.flow;
            flow_fn = fieldnames(flow);
            
            for i = 1:numel(visc_fn)
                
                field = visc_fn{i};
                
                if any(strcmp(field, flow_fn))
                    
                    try
                    obj.(field) = flow.(field);
                    catch
                    end
                end
            end
        end
    end
end

classdef Aerodynamics
    
    properties
        
        flightstate
        viscous
        stability
        object
        Aref = 1
        
        shield = false
        del
        impact
        shadow
        Cp
        Medge
        Vedge
        Pedge
        Tedge
        rho
        cf
        
        max_theta_beta_mach
        cone_table
        vmax
        
        impact_method = struct(...
            'newtonian', ["nose", "default"], ...
            'tangent_oblique_shock', ["foil", "tail", "wing", "wedge"], ...
            'tangent_cone', ["body", "fuse"], ...
            'oblique_shock_prandtl_meyer', ["cone"]);
        shadow_method = struct(...
            'oblique_shock_prandtl_meyer', ["fuse", "body", "foil", "wing"], ...
            'base', []); % "foil", "wing", "tail", "fuse", "body", "nose", "default"
    end
    
    methods
        function self = Aerodynamics(state, viscous, stability, Aref)
            
            if nargin > 0
                
                gamma = state.gamma;
                self.cone_table = tangent_cone_table(state.Minf, gamma);
                self.vmax = pi/2 * (((gamma + 1)/(gamma - 1))^0.5 - 1);
                self.flightstate = state;
                
                if nargin >= 2 && ~isempty(viscous)
                    
                    self.viscous = viscous;
                end
                if nargin >= 3 && ~isempty(stability)
                    
                    self.stability = stability;
                end
                if nargin >= 4 && ~isempty(Aref)
                    
                    self.Aref = Aref;
                end
            end
        end
        function [force_coeffs, forces, out] = trim(...
                self, weight, varargin)
            %% TODO: When integrating Structures module, this will have to be abstracted from Aerodynamics, so that it can be used for Aerostructual analysis also
            %TRIM attempts to find the trim angle of attack for an input
            %vehicle
            %   This is a simple routine that assumes no control surfaces,
            %   instead searching for the overall angle of attack at which 
            %   lift equals weight
            %
            %   Inputs:
            %   Aref - reference area to normalise force coefficients and
            %       calculate lifting force (m^2)
            %   weight - vehicle weight, must be in Newtons to compare to
            %       lifting force (N)
            %   varargin - objects that make up vehicle
            
            iter = 0;
            max_abs_alpha = deg2rad(45);
            
            % Set initial alpha to reasonable trim angle
            alpha = deg2rad(5);
            % Initial alpha step, applied positively if lift > weight,
            % otherwise negative
            initial_step = deg2rad(5);
            % Set initial alpha to reasonable trim angle
            self.flightstate.alpha = alpha;
            [force_coeffs, forces] = self.run(varargin{:});
            
            if isinf(weight)
                
                out.alpha_opt = inf;
                out.dCm_dalpha = inf;
                return
            end
            
            % Trimming loop
            while true
                
                lift = forces.L;
                LW_diff = lift - weight;
                
                if iter == 0
                    
                    alpha = alpha - (sign(LW_diff) * initial_step);
                else                
                    % Exit criteria:
                    %   Trim convergence (lift ~= weight)
                    %   Angle of attack is not longer changing due to trim 
                    %   or stuck at min/max alpha floor/ceiling
                    %   Cannot converge within tolerance
                    y = [prev_alpha, alpha];
                    
                    trimmed = abs(LW_diff/((lift + weight)/2)) < 0.01;
                    alpha_conv = ~diff(round(y, 6));
                    quit = (iter >= 50 && lift > weight);
                    
                    if trimmed || alpha_conv || quit
                        
                        % dCm = (Cm - previous.Cm)/(t.alpha - previous.alpha);
                        dCm = (force_coeffs.Cm - prev_force_co.Cm)/...
                            (force_coeffs.Cl - prev_force_co.Cl);
                        
                        if isnan(dCm), dCm = 1; end
                        
                        out.alpha_opt = self.flightstate.alpha;
                        out.dCm_dalpha = dCm;
                        break
                    else
                        %% TODO: Should use a simple line search here
                        % lift - weight gives a convex problem to be solved
                        % Currently assuming lift curve slope is rougly 
                        % linear, interpolating being most recent two runs 
                        % should always get a closer value
                        
                        % Interpolating to find alpha at which lift is
                        % equal to weight
                        x = [prev_lift, lift];
                        alpha = vec_interp(x(1), x(2), y(1), y(2), weight);
                        % Limiting output alpha within min/max bounds
                        alpha = max(...
                            min(alpha, max_abs_alpha), ...
                            -max_abs_alpha);
                    end
                end
                
                % Save previous values
                prev_force_co = force_coeffs;
                prev_lift = lift;
                prev_alpha = self.flightstate.alpha;
                
                % Update alpha and re-run
                self.flightstate.alpha = alpha;
                [force_coeffs, forces] = self.run(varargin{:});
                iter = iter + 1;
            end
        end
        function [dCm_dalpha, dCl_dbeta] = get_stability(...
                self, Cm, Cl, varargin)
            %GET_STABILITY calculates longitudinal and lateral stability
            %derivatives
            %   The current values of alpha and beta will be perturbed to
            %   produce the stability derivatives, which is useful when
            %   a trim routine has been carried out, providing stability
            %   derivatives for the trimmed state.
            %
            %   Inputs:
            %   Cm - the pitching moment coefficient for the current
            %       unperturbed flightstate
            %   Cl - the lift coeffient for the current unperturbed
            %       flightstate
            %   varargin - the input objects to be analysed
            %
            %   If only Cm is entered, only the longitudinal coefficient
            %   will be computed.
            %   If only Cl is entered, only the lateral coefficient will be
            %   computed.
            %   If neither Cm or Cl is entered, it is assumed that the
            %   current flightstate has not yet been assessed, and
            %   therefore will be assessed before calcuating stability
            %   derivatives.
            
            if isempty(Cm) && isempty(Cl)
                
                fc = self.run(varargin{:});
                Cm = fc.Cm;
                Cl = fc.Cl;
            end
            
            dCm_dalpha = [];
            dCl_dbeta = [];
            
            if ~isempty(Cm)
                %% Longitudinal stability
                % Setting up angle delta for angle of attack and sideslip
                angle_perturabtion = deg2rad(0.5);
                
                original_alpha = self.flightstate.alpha;
                perturbed_alpha = original_alpha + angle_perturabtion;
                self.flightstate.alpha = perturbed_alpha;
                
                pert_force_coeffs = self.run(varargin{:});
                
                dCm_dalpha = (pert_force_coeffs.Cm - Cm)/...
                    (perturbed_alpha - original_alpha);
                
                % Resetting alpha
                self.flightstate.alpha = original_alpha;
            end
            
            if ~isempty(Cl)
                %% Lateral stability
                original_beta = self.flightstate.beta;

                if original_beta == 0
                    % Multiplying by two to get full vehicle lift coefficient
                    Cl = Cl * 2;
                end

                perturbed_beta = original_beta + angle_perturabtion;
                self.flightstate.beta = perturbed_beta;
                
                %% TODO: Move this to subfunction which is called whenever beta ~= 0
                for i = numel(varargin):-1:1

                    part = varargin{i};

                    if ~isempty(part)

                        mirror = part;

                        mirror.x = flip(mirror.x, 2);
                        mirror.y = flip(-mirror.y, 2);
                        mirror.z = flip(mirror.z, 2);
                        mirror = mirror.get_data;

                        mirrored_part{i} = mirror;
                    end
                end

                fc1 = self.run(varargin{:});
                fc2 = self.run(mirrored_part{:});
                
                %%
                perturbed_Cl = fc1.Cl + fc2.Cl;

                % Resetting beta
                self.flightstate.beta = original_beta;

                dCl_dbeta = (perturbed_Cl - Cl)/...
                    (perturbed_beta - original_beta);
            end
        end
        function [outputs] = run(self, varargin)
            %RUN wraps the main function to allow multiple objects to be
            %analysed and aggregated, producing overall vehicle
            %characteristics
            for i = numel(varargin):-1:1
                
                part = varargin{i};
                part = part.get_data();
                
                [~, part_outputs(i)] = self.main(part);
            end
            
            if numel(varargin) == 1
                
                outputs = part_outputs;
                
                combined_fco = [outputs.fco];
                combined_forces = [outputs.forces];
                force_coeffs = combined_fco.sum();
                forces = combined_forces.sum();
            end
        end
        
        function [self, outputs] = main(self, part)
            
            self.object = part;
            
            state = self.flightstate;
            
            Unorm = state.Uvec/state.Uinf;
            
            data = part.data;
            unit_norm = data.unit_norm;
            area = data.area;
            
            dim = size(area);
            
            self = self.surface_velocity(unit_norm);
            d = asin(min(max(dotmat(-Unorm, unit_norm), -1), 1));
            d(area == 0) = 0;
            
            % Sanity check
            if any(isnan(d))
                
                error("inclination NaN");
            end
            
            self.del = d;
            
            %% Impact zone
            [im, sh, self.Cp, self.Medge, self.Pedge, self.Tedge, self.rho] = ...
                deal(zeros(dim));
            
            [~, ~, iid] = Aerodynamics.get_method(class(part), self.impact_method);
            % Panel faces flow as has non-zero area
            im(self.del > 0 & area > 0) = iid;
            
            % Angle to flow too large for attached shock
            % 99% used to avoid close calls being returned as NaN from
            % get_oblique_shock_angle
            im(im & self.del >= state.max_del * 0.99) = 1;
            
            % Any method other than newtonian requires attached shock
            shock = im > 1;
            target = tan(self.del(shock));
            
            % Get shockwave angles, any NaN values should have been dealt
            % with already, but setting to Newtonian just in case
            beta = self.get_oblique_shock_angle(target);
            im(shock(isnan(beta))) = 1;
            
            %% Shadow zone
            [~, ~, sid] = Aerodynamics.get_method(class(part), self.shadow_method);
            sh(self.del <= 0 & area > 0) = sid;
            
            % For non-conical parts, where the streamwise and
            % x-dimension are assumed to be parallel. First shaded panel
            % results in the rest of the strip being defined as shadow
            % region
            if ~part.conical
                
                for i = 1:dim(2)
                    
                    id = find(~im(:,i), 1, 'first');
                    if isempty(id), continue; end
                    
                    im(id:end, i) = 0;
                    sh(id:end, i) = sid;
                end
            end
            
            self.impact = im;
            self.shadow = sh;
            
            self = self.method_switcher(self.impact_method, im);
            
            % Mach_q essentially a placeholder
            % Replaced with interpolation from Mach = 0 normal to flow, to
            % first panel that does not use newtonian
            for i = 1:dim(2)
                
                int = self.Medge(:,i) == state.Mach_q | isnan(self.Medge(:,i));
                row = find(~int, 1);
                %% TODO: what if purely newtonian?
                x0 = self.del(row, i);
                y0 = self.Medge(row, i);
                
                if isempty(row) || y0 <= 0
                    x0 = 0;
                    y0 = self.flightstate.Minf;
                end
                
                x1 = pi/2;
                xq = self.del(int, i);
                y1 = 0;
                
                self.Medge(int, i) = y0 + (xq - x0).*((y1 - y0)./(x1 - x0));
            end
            
            self = self.method_switcher(self.shadow_method, sh);
            
            %% Run viscous analysis if viscous object supplied
            if ~isempty(self.viscous)
                % Copy necessary data from current object for analysis
                self.viscous = copy_class_properties(...
                    self, self.viscous);
                
                self.viscous = self.viscous.run();
                % Copying friction coefficient back to main object
                self.cf = self.viscous.cf;
            end
            
            outputs = self.get_outputs();
        end
        function [fco, forces] = get_force_data(self)
            
            fco = ForceCoeffs(...
                self.Cp, ...
                self.cf, ...
                self.object.data.area, ...
                self.object.data.unit_norm, ...
                self.flightstate.alpha, ...
                self.Aref);
            
            forces = Forces(fco, self.flightstate.q, self.Aref);
        end
        function outputs = get_outputs(self)
            %GET_OUTPUTS return a structure of chosen calculated parameters
            %   This avoids outputting the full object with lookup tables
            %   etc
            
            [fco, forces] = self.get_force_data();
            
            outputs = struct(...
                "force_coeffs", fco,...
                "forces", forces,...
                "Cp", self.Cp,...
                "cf", self.cf,...
                "viscous", self.viscous.get_outputs());
        end
        function self = method_switcher(self, methods, id)
            
            fn = fieldnames(methods);
            for k = 1:max(id(:))
                
                con = id == k;
                
                if any(con(:))
                    
                    self = self.(fn{k})(con);
                end
            end
        end
        function self = newtonian(self, con)
            
            state = self.flightstate;
            % theta = pi/2 - self.del(con);
            theta = self.del(con);
            
            Minf = state.Minf;
            Pinf = state.Pinf;
            g = state.gamma;
            
            a = 2./(g * Minf.^2);
            b = (((g + 1)^2) * Minf.^2)./((4 * g * Minf.^2) - (2 * (g - 1)));
            c = g/(g - 1);
            d = (1 - g + (2 * g * Minf.^2))./(g + 1);
            
            CpMax = a .* ((b.^c .* d) - 1);
            
            self.Cp(con) = CpMax .* sin(theta).^2;
            P = Pinf .* (1 + (self.Cp(con) * g * Minf.^2)/2);
            
            Vel = magmat(self.Vedge);
            
            Te = state.Tinf + (state.Uinf^2 - Vel(con).^2)/(2 * state.cp);
            
            self.rho(con) = P ./ (state.R * Te);
            self.Tedge(con) = Te;
            self.Pedge(con) = P;
            self.Medge(con) = state.Mach_q * ones(size(theta));
            
            % due_dx = state.q/(rhoe * ue) * cos(theta) * sin(theta);
        end
        
        function self = tangent_oblique_shock(self, con)
            
            theta = self.del(con);
            state = self.flightstate;
            
            Minf = state.Minf;
            Pinf = state.Pinf;
            g = state.gamma;
            
            beta = self.get_oblique_shock_angle(tan(theta), Minf, state.max_beta);
            
            % Fundamentals of Aerodynamics, Anderson 2001 & Development of an
            % Aerodynamics Code for the Optimisation of Hypersonic Vehicles Jazra
            % & Smart 2009
            Mn1 = Minf * sin(beta);
            
            numer = 1 + ((g - 1)/2) * Mn1.^2;
            denom = g * (Mn1.^2) - (g - 1)/2;
            
            Mn2 = (numer./denom).^0.5;
            
            % Anderson2006 p43
            K = Minf * theta;
            Pratio = 1 + (((g * (g + 1))/4) * K.^2) + (g * K.^2 .* (((g + 1)/4).^2 + (1./(K.^2))).^0.5);
            rratio = ((g + 1) * (Mn1.^2))./(2 + (g - 1) * (Mn1.^2));
            
            M2 = Mn2./(sin(beta - theta));
            
            self.Medge(con) = M2;
            % Anderson2006 p47
            % M1 = Minf as tangent assumes series of impinging shockwaves
            self.Cp(con) = (1/(0.5 * g * (Minf^2))) * (Pratio - 1);
            self.Pedge(con) = Pratio * Pinf;
            self.rho(con) = rratio * state.rinf;
            % Eq of state
            self.Tedge(con) = (Pratio./rratio) * state.Tinf;
        end
        function self = tangent_cone(self, con)
            
            theta = self.del(con);
            M1 = self.flightstate.Minf;
            g = self.flightstate.gamma;
            
            % Development of an Aerodynamics Code for the Optimisation of
            % Hypersonic Vehicles - gives M2 > M1
            % Mix between On hypersonic flow past unyawed cone & Approximate
            % solutions for supersonic flow over wedges and cones
            
            tau = asin(sin(theta) .* (((g + 1)/2) + (1./((M1 * sin(theta)).^2))).^0.5);
            
            % numer = (Minf^2) .* (cos(tau).^2) .* (1+2*(tau-theta));
            % denom = 1 + ((gamma-1)/2) * (Minf^2) * ((sin(tau).^2) - 2*((tau-theta).^2) .* (cos(theta).^2));
            % M = (numer./denom).^0.5;
            
            [~, M2] = halfspace(tau, self.cone_table);
            M2(M2 > M1) = M1;
            
            %% TODO: Using Kbeta instead of K?
            % Anderson2006 p143
            frac1 = ((g + 1) * ((M1 * sin(theta)).^2) + 2)./ ...
                ((g - 1) * ((M1 * sin(theta)).^2) + 2);
            
            frac2 = ((g + 1)/2) + 1./((M1 * sin(theta)).^2);
            
            Mn1 = M1 * sin(tau);
            
            K = M1 * theta;
            Pratio = 1 + (((g * (g + 1))/4) * K.^2) + (g * K.^2 .* (((g + 1)/4).^2 + (1./(K.^2))).^0.5);
            %% TODO: Same as tangent-wedge?
            % Anderson2006 p38
            rratio = ((g + 1) * (Mn1.^2))./(2 + (g - 1) * (Mn1.^2));
            
            self.Medge(con) = M2;
            self.Cp(con) = sin(theta.^2) .* (1 + (frac1 .* log(frac2)));
            self.Pedge(con) = Pratio * self.flightstate.Pinf;
            self.rho(con) = rratio * self.flightstate.rinf;
            % Eq of state
            self.Tedge(con) = Pratio./rratio .* self.flightstate.Tinf;
        end
        function self = oblique_shock_prandtl_meyer(self, con)
            
            state = self.flightstate;
            Minf = state.Minf;
            Pinf = state.Pinf;
            Tinf = state.Tinf;
            rinf = state.rinf;
            g = state.gamma;
            
            d = self.del;
            ddiff = [d(1,:); diff(d)];
            [row, col] = size(ddiff);
            
            %% TODO: Just have self.Value set inside functions, set V1 = self.V(i,:) in loop?
            Cp2 = self.Cp;
            M2 = self.Medge;
            P2 = self.Pedge;
            T2 = self.Tedge;
            r2 = self.rho;
            
            for i = 1:row
                
                coni = con(i,:);
                
                if ~any(coni), continue; end
                
                if i == 1
                    
                    M1 = repmat(Minf, 1, col);
                    P1 = repmat(Pinf, 1, col);
                    T1 = repmat(Tinf, 1, col);
                    r1 = repmat(rinf, 1, col);
                else
                    M1 = M2(i-1, :);
                    P1 = P2(i-1, :);
                    T1 = T2(i-1, :);
                    r1 = r2(i-1, :);
                end
                
                newt = ddiff(i,:) * 0.95 > state.max_del;
                shock = ddiff(i,:) > 1e-7 & ~newt;
                obs = coni & shock;
                
                if any(newt)
                    
                    col_arr = 0:col-1;
                    id_newt = i + (col_arr(newt))*row;
                    self = self.newtonian(id_newt);
                    Cp2(i, newt) = self.Cp(id_newt);
                    M2(i, newt) = self.Medge(id_newt);
                    P2(i, newt) = self.Pedge(id_newt);
                    T2(i, newt) = self.Tedge(id_newt);
                    r2(i, newt) = self.rho(id_newt);
                end
                
                if any(obs)
                    
                    [Cp2(i, obs), M2(i, obs), P2(i, obs), T2(i, obs), r2(i, obs)] = ...
                        oblique_shock(ddiff(i, obs), M1(obs), P1(obs), T1(obs), r1(obs));
                end
                
                %% Setting non-physical values to previous panel values
                set_prev = ~isfinite(M2(i,:)) | ~isreal(M2(i,:));
                
                if any(set_prev) && i > 1
                    
                    Cp2(i, set_prev) = Cp2(i-1, set_prev);
                    M2(i, set_prev) = M1(set_prev);
                    P2(i, set_prev) = P1(set_prev);
                    T2(i, set_prev) = T1(set_prev);
                    r2(i, set_prev) = r1(set_prev);
                end
                
                %% Prandtl-Meyer expansion
                % Second condition ensures shadow panels that have impact
                % panels further down the stream will be dealt with here
                % Otherwise impact panels would use uninitialised flow
                % properties
                % Must use shadow (~impact) instead of ~coni as using
                % latter could overwrite flow properties already calculated
                % from other impact methods
                pm = (coni & ~shock) | (self.shadow(i,:) & any(con(i+1:end, :)));
                
                if any(pm)
                    
                    [Cp2(i, pm), M2(i, pm), P2(i, pm), T2(i, pm), r2(i, pm)] = ...
                        prandtl_meyer(ddiff(i, pm), M1(pm), P1(pm), T1(pm), r1(pm));
                end
                
                %% TODO: Reattachment for large bases/reinclination to flow
                % Can't have one liner, sometimes alpha so high that
                % leading edge is set to base
                base = isnan(M2(i,:));
                
                if i > 1
                    
                    % If pm failed or if previous panel is base, next is also
                    base = base | Cp2(i-1,:) == self.get_base_pressure(Minf, g);
                end
                
                if any(base)
                    
                    col_arr = 0:col-1;
                    id_base = i + (col_arr(base))*row;
                    self = self.base(id_base);
                    Cp2(i, base) = self.Cp(id_base);
                    M2(i, base) = self.Medge(id_base);
                    P2(i, base) = self.Pedge(id_base);
                    T2(i, base) = self.Tedge(id_base);
                    r2(i, base) = self.rho(id_base);
                end
            end
            
            % Only return panel results asked for
            self.Cp = Cp2;
            self.Medge = M2;
            self.Pedge = P2;
            self.Tedge = T2;
            self.rho = r2;
            
            function [Cp2, M2, P2, T2, r2] = oblique_shock(del, M1, P1, T1, r1)
                
                del = del(:)';
                
                if nargin < 2 || isempty(M1), M1 = Minf; end
                if nargin < 3 || isempty(P1), P1 = Pinf; end
                if nargin < 4 || isempty(T1), T1 = Tinf; end
                if nargin < 5 || isempty(r1), r1 = rinf; end
                
                beta = self.get_oblique_shock_angle(tan(del), M1);
                
                M2 = (1./(sin(beta - del).^2) .* ...
                    ((g - 1) .* M1.^2 .* (sin(beta).^2) + 2)./ ...
                    (2 * g * (M1.^2) .* (sin(beta).^2) - (g - 1))).^0.5;
                
                % test = ((1 + ((g - 1)/2) * ((M1 .* sin(beta)).^2))./(g * ((M1 .* sin(beta)).^2) - (g - 1)/2)).^0.5;
                
                % Slight numerical error between two versions
                % Anderson2006 p38
                P2_P1 = 1 + ((2 * g)/(g + 1)) .* ...
                    ((M1.^2) .* (sin(beta).^2) - 1);
                % Where? Same as above?
                % P2 = P1 .* (2 * g * (M1.^2) .* ...
                % (sin(beta).^2) - (g - 1))/(g + 1);
                
                % Anderson2006 p38
                r2_r1 = ((g + 1) * (M1.^2) .* (sin(beta).^2))./ ...
                    (2 + (g - 1) * (M1.^2) .* (sin(beta).^2));
                
                P2 = P2_P1 .* P1;
                r2 = r2_r1 .* r1;
                % Eq of state
                T2 = P2_P1./r2_r1 .* T1;
                
                Cp2 = 4/(g + 1) * (sin(beta).^2) - 1./(M1.^2);
            end
            function [Cp2, M2, P2, T2, r2] = prandtl_meyer(dTheta, M1, P1, T1, r1)
                
                if nargin < 2 || isempty(M1), M1 = Minf; end
                if nargin < 3 || isempty(P1), P1 = Pinf; end
                if nargin < 4 || isempty(T1), T1 = Tinf; end
                if nargin < 5 || isempty(r1), r1 = rinf; end
                
                max_mach = 1000;
                M1(M1 > max_mach) = max_mach;
                
                pm_fun = @(M1) ((g + 1)/(g - 1)).^0.5 ...
                    * atan((((g - 1)/(g + 1)).*(M1.^2 - 1)).^0.5) ...
                    - atan(((M1.^2) - 1).^0.5);
                
                % Anderson p45
                vMp1 = pm_fun(M1);
                % Absolute value as del is relative to flow direction, not
                % convex/concavity
                vMp2 = vMp1 + abs(dTheta);
                
                [M2,~,~] = bisection(pm_fun, M1, 2 * max_mach, vMp2);
                
                [T2_T1, P2_P1, r2_r1] = self.isentropic(M1, M2);
                
                T2 = T2_T1 .* T1;
                P2 = P2_P1 .* P1;
                r2 = r2_r1 .* r1;
                
                Cp2 = 2 * ((P2/Pinf) - 1)/(g * (Minf^2));
            end
        end
        function beta = get_oblique_shock_angle(self, target, M1, ub)
            
            state = self.flightstate;
            
            if nargin < 3 || isempty(M1), M1 = state.Minf; end
            if nargin < 4 || isempty(ub)
                
                con = M1 == state.Minf;
                ub(con) = state.max_beta;
                
                if ~all(con)
                    
                    [~, ub(~con)] = halfspace(M1(~con), state.max_theta_beta_mach(:,[1 3]));
                end
            end
            
            g = state.gamma;
            
            % Has to be made here due to M1 dependency
            oblque_shock_eqn = @(B) 2*cot(B) .* ((M1.^2) .* (sin(B).^2) - 1)./((M1.^2) .* (g + cos(2*B)) + 2);
            [beta, ~, ~] = bisection(oblque_shock_eqn, 0, ub, target, 1e-3);
        end
        function self = base(self, con)
            
            state = self.flightstate;
            %             self.Medge(con) = state.Minf;
            %             self.Tedge(con) = state.Tinf;
            %             self.Pedge(con) = state.Pinf;
            %             self.rho(con) = state.rinf;
            
            self.Medge(con) = 0;
            self.Tedge(con) = 0;
            self.Pedge(con) = 0;
            self.rho(con) = 0;
            self.Cp(con) = self.get_base_pressure(state.Minf, state.gamma);
        end
        function [self, tang, surf] = surface_velocity(self, unorm, normalise)
            
            tang = crossmat(unorm, self.flightstate.Uvec);
            surf = crossmat(tang, unorm);
            
            self.Vedge = surf;
            
            if nargin >= 3 && normalise
                
                tang_mag = magmat(tang);
                surf_mag = magmat(surf);
                tang = tang./tang_mag;
                surf = surf./surf_mag;
            end
        end
        function [T2_T1, P2_P1, r2_r1] = isentropic(self, M1, M2)
            % Isentropic flow relations
            %% TODO: Move to flightstate?
            g = self.flightstate.gamma;
            
            T2_T1 = (1 + ((g - 1)/2) * (M1.^2)) ...
                ./ (1 + ((g - 1)/2) * (M2.^2));
            
            P2_P1 = T2_T1.^(g/(g - 1));
            r2_r1 = T2_T1.^(1/(g - 1));
        end
        function centre_of_pressure = get_centre_of_pressure(self, centre)
            
            xyz_Cp = reshape(centre .* self.Cp, [], 3, 1);
            sum_Cp = sum(self.Cp(:));
            
            centre_of_pressure = xyz_Cp/sum_Cp;
        end
    end
    
    methods (Static)
        function [fname, res, i] = get_method(name, self)
            
            fields = fieldnames(self);
            
            % In reverse order so that default is first method
            for i = numel(fields):-1:1
                
                fname = fields{i};
                res = self.(fname);
                
                if isstring(res) && any(contains(name, res, 'IgnoreCase', true))
                    return
                end
            end
        end
        function Cp = get_base_pressure(M, g)
            
            if nargin == 1
                
                Cp = -1./(M.^2);
            else
                Cp = (2./(g * M.^2)) .* ((2/(g + 1))^1.4 .* (1./M).^2.8 .* ...
                    ((2 * g * (M.^2) - (g - 1))/(g + 1)) - 1);
            end
        end
    end
end
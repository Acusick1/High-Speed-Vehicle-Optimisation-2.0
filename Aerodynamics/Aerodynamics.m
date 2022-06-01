classdef Aerodynamics
    
    properties
        
        flow
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
        force_coeffs
        forces
        visc
        CoP
        maxTable
        coneTable
        vmax
        dCm_dalpha
        dCl_dbeta
        
        trim = false
        part
        weight
        alpha_opt
        latStab = true
    end
    
    properties (Dependent)
        
    end
    
    properties %(Constant)
        
        impact_method = struct(...
            'newtonian', ["nose", "default"], ...
            'tangentobs', ["foil", "tail", "wing", "wedge"], ...
            'tangentcone', ["body", "fuse"], ...
            'obspm', ["cone"]);
        shadow_method = struct(...
            'obspm', ["fuse", "body", "foil", "wing"], ...
            'base', []); % "foil", "wing", "tail", "fuse", "body", "nose", "default"
    end
    
    methods
        function obj = Aerodynamics(viscous, trim, states, nParts)
            
            if nargin > 0
                
                temp = Aerodynamics();
                
                if nargin > 1 && viscous, temp.visc = Viscous(); end
                if nargin < 2 || isempty(trim), trim = false; end
                
                temp.trim = trim;
                [~, temp.maxTable] = theta_beta_mach_curves();
                
                if nargin >= 3
                    
                    for i = numel(states):-1:1

                        obj(i,:) = temp;
                        obj(i,:) = obj(i,:).set_flow(states(i));
                    end
                    if nargin >= 4 && ~isempty(nParts)
                        
                        obj = repelem(obj, 1, nParts);
                    end
                else
                    obj = temp;
                end
            end
        end
        function obj = set_flow(obj, state)
            
            g = state.gamma;
            obj.coneTable = tangent_cone_table(state.Minf, g);
            obj.vmax = pi/2 * (((g + 1)/(g - 1))^0.5 - 1);
            obj.flow = state;
        end
        function obj = analyse(obj, varargin)
            
            [nStates, nParts] = size(obj);
            
            for i = 1:nStates
                
                temp = obj(i,1);
                W = temp.weight;
                
                % If next state has same Mach and altitude, do not trim.
                % Only reason this occurs is if uncertainty is desired in
                % trim angle. Therefore do not re-trim, just apply 
                % uncertainty 
                if i > 1 && temp.trim && temp.flow.Minf == obj(i-1,1).flow.Minf
                    
                    temp.trim = false;
                    
                    for j = 1:nParts

                        % If vehicle unable to trim, set alpha equal to
                        % previous for size consistency
                        try
                            obj(i,j).flow.alpha = t.alpha + ...
                                temp.flow.alpha;
                        catch
                            obj(i,j).flow.alpha = obj(i-1,j).flow.alpha + ...
                                temp.flow.alpha;
                        end
                        
                        obj(i,j).alpha_opt = obj(i-1,j).alpha_opt;
                        obj(i,j).dCm_dalpha = obj(i-1,j).dCm_dalpha;
                        obj(i,j).dCl_dbeta = obj(i-1,j).dCl_dbeta;
                    end
                end
                
                t.count = 1;
                p.Cm = nan;
                p.alpha = nan;
                alpha_conv = false;
                
                while true
                    
                    for j = 1:nParts
                        
                        geometry = varargin{j};
                        
                        if ~isempty(geometry)
                            
                            obj(i,j) = obj(i,j).main(geometry);
                        end
                    end
                    
                    a = obj(i,1).flow.alpha;
                    f = [obj(i,:).forces];
                    fc = [obj(i,:).force_coeffs];
                    lift = sum([f.L]);
                    Cm = sum([fc.Cm]);
                    Cl = sum([fc.Cl]);
                    
                    if ~temp.trim
                        break
                    else
                        
                        % Trim convergence criteria:
                        % Actual convergence ie. lift ~= weight,
                        % Cannot converge within tolerance
                        trimmed = 0.99 * W < lift && ...
                            1.01 * W > lift;
                        
                        quit = t.count >= 50 && lift > W;
                        
                        if t.count > 1 && (trimmed || alpha_conv || quit)
                            
                            if t.alpha > p.alpha
                                
                                % dCm = (Cm - p.Cm)/(t.alpha - p.alpha);
                                dCm = (Cm - p.Cm)/(Cl - p.Cl);
                            else
                                % dCm = (p.Cm - Cm)/(p.alpha - t.alpha);
                                dCm = (p.Cm - Cm)/(p.Cl - Cl);
                            end
                            
                            if isnan(dCm), dCm = 1; end
                            obj(i,1).dCm_dalpha = dCm;
                            
                            break
                        else
                            p.Cm = Cm;
                            p.Cl = Cl;
                            p.alpha = obj(i,1).flow.alpha;
                            
                            [t, alpha_conv] = update_alpha(lift, t);
                        end
                        
                        %% Next loop setup
                        for j = 1:nParts
                            
                            obj(i,j).flow.alpha = t.alpha;
                            obj(i,j).alpha_opt = t.alpha;
                        end
                    end
                end
                
                %% TODO: Better way of only coming in here for non-uncertainty flight states
                if temp.latStab && temp.trim
                    
                    for j = nParts:-1:1
                        
                        side1(i,j) = obj(i,j);
                        side2(i,j) = obj(i,j);
                        side1(i,j).flow.beta = deg2rad(0.5);
                        side2(i,j).flow.beta = deg2rad(0.5);
                        
                        geometry = varargin{j};
                        
                        if ~isempty(geometry)
                            
                            mirror = geometry;
                            
                            mirror.x = flip(mirror.x, 2);
                            mirror.y = flip(-mirror.y, 2);
                            mirror.z = flip(mirror.z, 2);
                            mirror = mirror.get_data;
                            
                            side1(i,j) = side1(i,j).main(geometry);
                            side2(i,j) = side2(i,j).main(mirror);
                        end
                    end
                    
                    fc1 = [side1(i,:).force_coeffs];
                    fc2 = [side2(i,:).force_coeffs];
                    Clb = sum([fc1.Cl] + [fc2.Cl])/2;
                    
                    obj(i,1).dCl_dbeta = (Clb - Cl)/deg2rad(0.5);
                end
            end
            
            function [t, conv] = update_alpha(lift, t)
                
                LWdiff = lift - W;
                
                if t.count == 1
                    
                    % Can be temp here as only coming in initial instance
                    t.alpha = temp.flow.alpha;
                    
                    if LWdiff > 0
                        
                        x = [nan LWdiff];
                        y = [nan t.alpha];
                        alpha = t.alpha - deg2rad(5);
                    else
                        x = [LWdiff nan];
                        y = [t.alpha nan];
                        alpha = t.alpha + deg2rad(5);
                    end
                else
                    alpha = t.alpha;
                    x = t.x;
                    y = t.y;
                    
                    if t.count == 2
                        
                        x(isnan(x)) = LWdiff;
                        y(isnan(y)) = alpha;
                    else
                        %% TODO: Cleaner logic
                        if (LWdiff < 0 && LWdiff > x(1)) || all(x(1) > [0 LWdiff])
                            % If new < 0 and > lower bound, or
                            % lower bound > 0 and new < lower bound
                            % lower bound = new
                            x(1) = LWdiff;
                            y(1) = alpha;
                              
                        elseif LWdiff > 0 && LWdiff < x(2) || all(x(2) < [0 LWdiff])
                            % If new > 0 and < upper bound, or
                            % upper bound < 0 and new > upper bound
                            % upper bound = new
                            x(2) = LWdiff;
                            y(2) = alpha;
                        end
                    end
                    
                    alpha = y(1) - (x(1) * (y(2) - y(1))/(x(2) - x(1)));
                    alpha = max(min(alpha, deg2rad(45)), deg2rad(-45));
                end
                
                % Comparing current and next alphas
                alphas = [t.alpha, alpha];
                
                % If prev and next angle of attack are too high
                too_steep = all(abs(alphas) > deg2rad(35));
                % If alpha is barely changing
                approx_conv = ~diff(round(alphas, 6));
                
                conv = too_steep || approx_conv;
                
                if ~conv
                    
                    t.alpha = alpha;
                    t.x = x;
                    t.y = y;
                end
                
                t.count = t.count + 1;
            end
        end
        function obj = main(obj, part)
            
            f = obj.flow;
            
            Unorm = f.U/f.Uinf;
            
            if isa(part, 'Aerofoil')
                
                data = part.data;
                unit_norm = data.unit_norm;
                area = data.area;
            else
                data = part.quad_data;
                unit_norm = data.norm./data.mag;
                area = data.area;
            end
            
            dim = size(area);
            
            obj = obj.surface_velocity(unit_norm);
            d = asin(min(max(dotmat(-Unorm, unit_norm), -1), 1));
            d(area == 0) = 0;
            
            con = isnan(d);
            if any(con(:)), error("inclination NaN"); end
            
            d(con) = 0;
            obj.del = d;
            
            %% Impact zone
            [im, sh, obj.Cp, obj.Medge, obj.Pedge, obj.Tedge, obj.rho] = deal(zeros(dim));
            
            [~, ~, iid] = Aerodynamics.get_method(class(part), obj.impact_method);
            % Panel faces flow as has non-zero area
            im(obj.del > 0 & area > 0) = iid;
            % im(obj.del > 0) = iid;
            % Angle to flow too large for attached shock
            % 99% used to avoid close calls being returned as NaN from tbm
            im(im & obj.del >= f.maxDel * 0.99) = 1;
            
            % Any method other than newtonian requires attached shock
            shock = im > 1;
            target = tan(obj.del(shock));
            
            beta = obj.tbm(target);
            rm = isnan(beta);
            
            %% Potential error
            % If this ever shows up, obj.del > maxDel is not enough to
            % flag all detached shockwaves. Add below back in:
            im(shock(rm)) = 1;
            if any(rm)
                
                error("FIX BETA NAN")
            end
            
            %% Shadow zone
            [~, ~, sid] = Aerodynamics.get_method(class(part), obj.shadow_method);
            sh(obj.del <= 0 & area > 0) = sid;

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
            
            obj.impact = im;
            obj.shadow = sh;

            obj = obj.method_switcher(obj.impact_method, im);
            
            % Machq essentially a placeholder
            % Replaced with interpolation from Mach = 0 normal to flow, to
            % first panel that does not use newtonian
            for i = 1:dim(2)
                
                int = obj.Medge(:,i) == f.Machq | isnan(obj.Medge(:,i));
                a = find(~int, 1);
                %% TODO: what if purely newtonian?
                x0 = obj.del(a, i);
                y0 = obj.Medge(a, i);
                
                if isempty(a) || y0 <= 0    
                    x0 = 0;
                    y0 = obj.flow.Minf;
                end
              
                x1 = pi/2;
                xq = obj.del(int, i);
                y1 = 0;
                
                obj.Medge(int, i) = y0 + (xq - x0).*((y1 - y0)./(x1 - x0));
            end
            
            obj = obj.method_switcher(obj.shadow_method, sh);
            
            obj.part = part;
            
            if ~isempty(obj.visc)
                
                v = Viscous.from_aerodynamics(obj, obj.visc);
                v = v.test();
                cf = v.cf;
                obj.visc = v;
            else
                cf = [];
            end
            
            obj.force_coeffs = ForceCoeffs(data, obj.Cp, cf, f.alpha, obj.Aref);
            obj.forces = Forces(obj.force_coeffs, f, obj.Aref);
        end
        function obj = method_switcher(obj, methods, id)
            
            fn = fieldnames(methods);
            for k = 1:max(id(:))
                
                con = id == k;
                
                if any(con(:))
                    
                    obj = obj.(fn{k})(con);
                end
            end
        end
        function a = get_centre_of_pressure(obj, centre)
            
            xyzCp = reshape(centre .* obj.Cp, [], 3, 1);
            sumCp = sum(obj.Cp(:));
            
            a = xyzCp/sumCp;
        end
        function obj = newtonian(obj, con)
            
            f = obj.flow;
            % theta = pi/2 - obj.del(con);
            theta = obj.del(con);
            
            Minf = f.Minf;
            Pinf = f.Pinf;
            g = f.gamma;
            
            a = 2./(g * Minf.^2);
            b = (((g + 1)^2) * Minf.^2)./((4 * g * Minf.^2) - (2 * (g - 1)));
            c = g/(g - 1);
            d = (1 - g + (2 * g * Minf.^2))./(g + 1);
            
            CpMax = a .* ((b.^c .* d) - 1);
            
            obj.Cp(con) = CpMax .* sin(theta).^2;
            P = Pinf .* (1 + (obj.Cp(con) * g * Minf.^2)/2);
            
            Vel = magmat(obj.Vedge);
            
            Te = f.Tinf + (f.Uinf^2 - Vel(con).^2)/(2 * f.cp);
            
            obj.rho(con) = P ./ (f.R * Te);
            obj.Tedge(con) = Te;
            obj.Pedge(con) = P;
            obj.Medge(con) = f.Machq * ones(size(theta));
            
            % due_dx = f.q/(rhoe * ue) * cos(theta) * sin(theta);
        end
        
        function obj = tangentobs(obj, con)
            
            theta = obj.del(con);
            f = obj.flow;
            
            Minf = f.Minf;
            Pinf = f.Pinf;
            g = f.gamma;
            
            beta = obj.tbm(tan(theta), Minf, f.maxBeta);
            
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
            
            obj.Medge(con) = M2;
            % Anderson2006 p47
            % M1 = Minf as tangent assumes series of impinging shockwaves
            obj.Cp(con) = (1/(0.5 * g * (Minf^2))) * (Pratio - 1);
            obj.Pedge(con) = Pratio * Pinf;
            obj.rho(con) = rratio * f.rinf;
            % Eq of state
            obj.Tedge(con) = (Pratio./rratio) * f.Tinf;
        end
        function obj = tangentcone(obj, con)

            theta = obj.del(con);
            M1 = obj.flow.Minf;
            g = obj.flow.gamma;
            
            % Development of an Aerodynamics Code for the Optimisation of
            % Hypersonic Vehicles - gives M2 > M1
            % Mix between On hypersonic flow past unyawed cone & Approximate
            % solutions for supersonic flow over wedges and cones
            
            tau = asin(sin(theta) .* (((g + 1)/2) + (1./((M1 * sin(theta)).^2))).^0.5);
            
            % numer = (Minf^2) .* (cos(tau).^2) .* (1+2*(tau-theta));
            % denom = 1 + ((gamma-1)/2) * (Minf^2) * ((sin(tau).^2) - 2*((tau-theta).^2) .* (cos(theta).^2));
            % M = (numer./denom).^0.5;
            
            [~, M] = halfspace(tau, obj.coneTable);
            M(M > M1) = M1;
            
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
            
            obj.Medge(con) = M;
            obj.Cp(con) = sin(theta.^2) .* (1 + (frac1 .* log(frac2)));
            obj.Pedge(con) = Pratio * obj.flow.Pinf;
            obj.rho(con) = rratio * obj.flow.rinf;
            % Eq of state
            obj.Tedge(con) = Pratio./rratio .* obj.flow.Tinf;
        end
        function obj = obspm(obj, con)
            
            f = obj.flow;
            Minf = f.Minf;
            Pinf = f.Pinf;
            Tinf = f.Tinf;
            rinf = f.rinf;
            g = f.gamma;
            
            d = obj.del;
            ddiff = [d(1,:); diff(d)];
            [row, col] = size(ddiff);
            
            %% TODO: Just have obj.Value set inside functions, set V1 = obj.V(i,:) in loop?
            Cp2 = obj.Cp;
            M2 = obj.Medge;
            P2 = obj.Pedge;
            T2 = obj.Tedge;
            r2 = obj.rho;
            
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
                
                newt = ddiff(i,:) * 0.95 > f.maxDel;
                shock = ddiff(i,:) > 1e-7 & ~newt;
                obs = coni & shock;
                
                if any(newt)
                    
                    col_arr = 0:col-1;
                    id_newt = i + (col_arr(newt))*row;
                    obj = obj.newtonian(id_newt);
                    Cp2(i, newt) = obj.Cp(id_newt);
                    M2(i, newt) = obj.Medge(id_newt);
                    P2(i, newt) = obj.Pedge(id_newt);
                    T2(i, newt) = obj.Tedge(id_newt);
                    r2(i, newt) = obj.rho(id_newt);
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
                pm = (coni & ~shock) | (obj.shadow(i,:) & any(con(i+1:end, :)));
                
                if any(pm)
                    
                    [Cp2(i, pm), M2(i, pm), P2(i, pm), T2(i, pm), r2(i, pm)] = ...
                        prandtlmeyer(ddiff(i, pm), M1(pm), P1(pm), T1(pm), r1(pm));
                end
                
                %% TODO: Reattachment for large bases/reinclination to flow
                % Can't have one liner, sometimes alpha so high that
                % leading edge is set to base
                base = isnan(M2(i,:));
                
                if i > 1
                    
                    % If pm failed or if previous panel is base, next is also
                    base = base | Cp2(i-1,:) == obj.base_pressure(Minf, g);
                end
                
                if any(base)
                    
                    col_arr = 0:col-1;
                    id_base = i + (col_arr(base))*row;
                    obj = obj.base(id_base);
                    Cp2(i, base) = obj.Cp(id_base);
                    M2(i, base) = obj.Medge(id_base);
                    P2(i, base) = obj.Pedge(id_base);
                    T2(i, base) = obj.Tedge(id_base);
                    r2(i, base) = obj.rho(id_base);
                end
            end
            
            % Only return panel results asked for
            obj.Cp = Cp2;
            obj.Medge = M2;
            obj.Pedge = P2;
            obj.Tedge = T2;
            obj.rho = r2;
            
            function [Cp2, M2, P2, T2, r2] = oblique_shock(del, M1, P1, T1, r1)
                
                del = del(:)';
                
                if nargin < 2 || isempty(M1), M1 = Minf; end
                if nargin < 3 || isempty(P1), P1 = Pinf; end
                if nargin < 4 || isempty(T1), T1 = Tinf; end
                if nargin < 5 || isempty(r1), r1 = rinf; end
                
                beta = obj.tbm(tan(del), M1);
                
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
            function [Cp2, M2, P2, T2, r2] = prandtlmeyer(dTheta, M1, P1, T1, r1)
                
                if nargin < 2 || isempty(M1), M1 = Minf; end
                if nargin < 3 || isempty(P1), P1 = Pinf; end
                if nargin < 4 || isempty(T1), T1 = Tinf; end
                if nargin < 5 || isempty(r1), r1 = rinf; end
                
                maxMach = 1000;
                M1(M1 > maxMach) = maxMach;
                
                pm_fun = @(M1) ((g + 1)/(g - 1)).^0.5 ...
                    * atan((((g - 1)/(g + 1)).*(M1.^2 - 1)).^0.5) ...
                    - atan(((M1.^2) - 1).^0.5);
                
                % Anderson p45
                vMp1 = pm_fun(M1);
                % Absolute value as del is relative to flow direction, not
                % convex/concavity
                vMp2 = vMp1 + abs(dTheta);
                
                [M2,~,~] = bisection(pm_fun, M1, 2*maxMach, vMp2);
                
                [T2_T1, P2_P1, r2_r1] = obj.isentropic(M1, M2);
                
                T2 = T2_T1 .* T1;
                P2 = P2_P1 .* P1;
                r2 = r2_r1 .* r1;
                
                Cp2 = 2 * ((P2/Pinf) - 1)/(g * (Minf^2));
            end
        end
        function beta = tbm(obj, target, M1, ub)
            
            f = obj.flow;
            
            if nargin < 3 || isempty(M1), M1 = f.Minf; end
            if nargin < 4 || isempty(ub)
                
                con = M1 == f.Minf;
                ub(con) = f.maxBeta;
                
                if ~all(con)
                    
                    [~, ub(~con)] = halfspace(M1(~con), obj.maxTable(:,[1 3]));
                end
            end
            
            g = f.gamma;
            
            % Has to be made here due to M1 dependency
            tbmFun = @(B) 2*cot(B) .* ((M1.^2) .* (sin(B).^2) - 1)./((M1.^2) .* (g + cos(2*B)) + 2);
            [beta, ~, ~] = bisection(tbmFun, 0, ub, target, 1e-3);
        end
        function obj = base(obj, con)

            f = obj.flow;
%             obj.Medge(con) = f.Minf;
%             obj.Tedge(con) = f.Tinf;
%             obj.Pedge(con) = f.Pinf;
%             obj.rho(con) = f.rinf;
            
            obj.Medge(con) = 0;
            obj.Tedge(con) = 0;
            obj.Pedge(con) = 0;
            obj.rho(con) = 0;
            obj.Cp(con) = obj.base_pressure(f.Minf, f.gamma);
        end
        function [obj, tang, surf] = surface_velocity(obj, unorm, normalise)
            
            tang = crossmat(unorm, obj.flow.U);
            surf = crossmat(tang, unorm);
            
            obj.Vedge = surf;
            
            if nargin >= 3 && normalise
                
                tang_mag = magmat(tang);
                surf_mag = magmat(surf);
                tang = tang./tang_mag;
                surf = surf./surf_mag;
            end
        end
        function [T2_T1, P2_P1, r2_r1] = isentropic(obj, M1, M2)
            %% Isentropic flow relations
            % Move to flightstate?
            g = obj.flow.gamma;
            
            T2_T1 = (1 + ((g - 1)/2) * (M1.^2)) ...
                ./ (1 + ((g - 1)/2) * (M2.^2));
            
            P2_P1 = T2_T1.^(g/(g - 1));
            r2_r1 = T2_T1.^(1/(g - 1));
        end
    end
    
    methods (Static)
        
        function [fname, res, i] = get_method(name, obj)
            
            fields = fieldnames(obj);
            
            % In reverse order so that default is first method
            for i = numel(fields):-1:1
                
                fname = fields{i};
                res = obj.(fname);
                
                if isstring(res) && any(contains(name, res, 'IgnoreCase', true))
                    return
                end
            end
        end
        function Cp = base_pressure(M, g)
            
            if nargin == 1
                
                Cp = -1./(M.^2);
            else
                Cp = (2./(g * M.^2)) .* ((2/(g + 1))^1.4 .* (1./M).^2.8 .* ...
                    ((2 * g * (M.^2) - (g - 1))/(g + 1)) - 1);
            end
        end
    end
end
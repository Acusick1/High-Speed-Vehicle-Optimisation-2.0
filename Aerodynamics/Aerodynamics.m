classdef Aerodynamics
    
    properties
        
        flow
        Aref
        shield = false
        viscous = false
        del
        impact
        shadow
        Cp
        Cf
        Mach
        V
        Pre
        T
        rho
        force_coeffs
        forces
        CoP
        maxTable
        pm_fun
        vmax
    end
    
    properties (Constant)
        
        impact_method = struct(...
            'newtonian', ["nose", "default", "body", "fuse", "wing"], ...
            'tangentobs', ["foil", "tail"], ...
            'tangentcone', [], ...
            'obspm', []);
        shadow_method = struct(...
            'obspm', [], ...
            'base', []); % "foil", "wing", "tail", "body", "nose", "default"
    end
    
    properties (Dependent)
        
    end
    
    methods
        function obj = Aerodynamics(viscous)
            
            if nargin > 0
                
                if nargin < 1 || isempty(viscous), viscous = false; end
                
                obj.viscous = viscous;
                [~, obj.maxTable] = theta_beta_mach_curves();
            end
        end
        function new_obj = analyse(obj, states, Aref, varargin)
            
            if isempty(Aref), Aref = 1; end
            
            for i = numel(states):-1:1
                %% TODO: Clean this hack up
                obj(i).maxTable = obj(1).maxTable;
                obj(i).flow = states(i);
                
                g = states(i).gamma;
                obj(i).pm_fun = @(M1) ((g + 1)/(g - 1)).^0.5 * atan((((g - 1)/(g + 1)).*(M1.^2 - 1)).^0.5) - atan(((M1.^2) - 1).^0.5);
                obj(i).Aref = Aref;
                obj(i).vmax = pi/2 * (((g + 1)/(g - 1))^0.5 - 1);
            end
            
            nStates = numel(obj);
            nParts = numel(varargin);
            
            new_obj = repelem(obj', 1, nParts);
            
            for i = 1:nStates
                
                for j = 1:nParts
                    
                    if ~isempty(varargin{j})
                        
                        new_obj(i,j) = new_obj(i,j).main(varargin{j});
                    end
                end
            end
        end
        function obj = main(obj, part)
            
            f = obj.flow;
            
            Unorm = f.U/f.Uinf;
            
            if isa(part, 'Aerofoil')
                %% TODO: Put in aerofoil.data
                upNorm = Aerofoil.normal([part.xu part.zu]);
                loNorm = -Aerofoil.normal([part.xl part.zl]);
                
                norm(:,:,1) = [upNorm(:,1) loNorm(:,1)];
                norm(:,:,3) = [upNorm(:,2) loNorm(:,2)];
                
                mag = magmat(norm);
                unit_norm = norm./mag;
                points = reshape([part.xu part.xl part.zu part.zl], [], 2, 2);
                area = magmat(diff(points, 1, 1));
                centre = (points(1:end-1,:,:) + points(2:end,:,:))/2;
                
                data.norm = norm;
                data.mag = mag;
                data.unitNorm = unit_norm;
                data.area = area;
                data.centre = centre;
            else
                data = part.quad_data;
                unit_norm = data.norm./data.mag;
                area = data.area;
                centre = data.centre;
            end
            
            dim = size(area);
            
            obj = obj.surface_velocity(unit_norm);
            d = asin(dotmat(-Unorm, unit_norm));
            d(area == 0) = 0;
            
            con = isnan(d);
            
            %% TODO: Does area == 0 remove need for below?
            
            if any(con)
                
                error("fix nan del")
            end
            
            obj.del = d;
            
            %% Impact zone
            
            [im, sh, obj.Cp, obj.Mach, obj.Pre, obj.T, obj.rho] = deal(zeros(dim));
            
            [~, ~, iid] = Aerodynamics.get_method(class(part), obj.impact_method);
            % Panel faces flow as has non-zero area
            im(obj.del > 0 & area > 0) = iid;
            % im(obj.del > 0) = iid;
            % Angle to flow too large for attached shock
            % 95% used to avoid close calls being returned as NaN from tbm
            im(im & obj.del >= f.maxDel * 0.95) = 1;
            
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
            
            %%
            [~, ~, sid] = Aerodynamics.get_method(class(part), obj.shadow_method);
            sh(obj.del <= 0 & area > 0) = sid;
            % sh(obj.del <= 0) = sid;
            
            obj.impact = im;
            obj.shadow = sh;
            
            %% TODO: Check these give the same answer, delete looper
            % obj = obj.looper();
            obj = obj.method_switcher(obj.impact_method, im);
            obj = obj.method_switcher(obj.shadow_method, sh);
            %%
            
            if obj.viscous
                
                obj.Cf = obj.eckert(centre, area, part.conical);
            end
            
            obj.force_coeffs = ForceCoeffs(data, obj.Cp, obj.Cf, f.alpha, obj.Aref);
            obj.forces = Forces(obj.force_coeffs, f);
        end
        function obj = looper(obj)
            % If decide to use this, have methods return values rather than
            % objects, avoid the con nonsense
            [row, col] = size(obj.del);
            con = zeros(row, col);
            
            im = obj.impact;
            sh = obj.shadow;
            
            for i = 1:row
                
                temp = con;
                
                if i > 1
                    % If not at leading edge, previous panel shadow method 
                    % non-Newtonian and Mach has not been set, set to base
                    base = sh(i,:) & obj.Mach(i-1, :) == 0;
                    im(i, base) = 0;
                    sh(i, base) = 2;
                end
                
                if any(im(i,:))
                    
                    temp(i,:) = im(i,:);
                    obj = obj.method_switcher(obj.impact_method, temp);
                end
                
                if any(sh(i,:))
                    
                    temp(i,:) = sh(i,:);
                    obj = obj.method_switcher(obj.shadow_method, temp);
                end
                
                if any(~isreal(obj.Cp) | ~isfinite(obj.Cp(:)))
                    
                    return
                end
            end
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
            
            Vel = magmat(obj.V);
            
            Te = f.Tinf + (f.Uinf^2 - Vel(con).^2)/(2 * f.cp);
            
            obj.rho(con) = P ./ (f.R * Te);
            obj.T(con) = Te;
            obj.Pre(con) = P;
            obj.Mach(con) = f.Machq * ones(size(theta));
            
            
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
            Pratio = 1 + ((2 * g)/(g + 1)) * ((Mn1.^2) - 1);
            
            M2 = Mn2./(sin(beta - theta));
            
            obj.Mach(con) = M2;
            % Assumed equation (Isentropic < flow is not)
            obj.Cp(con) = (1/(0.5 * g * (Minf^2))) * (Pratio - 1);
            % M1 = Minf as tangent assumes series of impinging shockwaves
            obj.Cp(con) = 4/(g + 1) * (sin(beta).^2) - 1./(Minf.^2);
            obj.Pre(con) = Pratio * Pinf;
        end
        function obj = tangentcone(obj, con)
            
            theta = obj.del(con);
            M1 = obj.flow.Minf;
            
            % Development of an Aerodynamics Code for the Optimisation of
            % Hypersonic Vehicles - gives M2 > M1
            % Mix between On hypersonic flow past unyawed cone & Approximate
            % solutions for supersonic flow over wedges and cones
            
            tau = asin(sin(theta) .* (((g + 1)/2) + (1./((M1 * sin(theta)).^2))).^0.5);
            
            % numer = (Minf^2) .* (cos(tau).^2) .* (1+2*(tau-theta));
            % denom = 1 + ((gamma-1)/2) * (Minf^2) * ((sin(tau).^2) - 2*((tau-theta).^2) .* (cos(theta).^2));
            % M = (numer./denom).^0.5;
            
            % Exact Taylor-Maccoll solution?
            for i = length(tau):-1:1
                
                [~,M,~] = solvecone(tau, M1, g);
            end
            
            frac1 = ((g + 1) * ((M1 * sin(theta)).^2) + 2)./ ...
                ((g - 1) * ((M1 * sin(theta)).^2) + 2);
            
            frac2 = ((g + 1)/2) + 1./((M1 * sin(theta)).^2);
            
            Pratio = (2 * g * M1.^2 * (sin(tau).^2) - (g - 1))/(g + 1);
            
            obj.Mach(con) = M;
            obj.Cp(con) = (theta.^2) .* (1 + (frac1 .* log(frac2)));
            obj.Pre(con) = Pratio * obj.flow.Pinf;
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
            
            Cp2 = obj.Cp;
            M2 = obj.Mach;
            P2 = obj.Pre;
            T2 = obj.T;
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
                
                %% TODO: newtonian
                
                if any(obs)
                    
                    [Cp2(i, obs), M2(i, obs), P2(i, obs), T2(i, obs), r2(i, obs)] = ...
                        oblique_shock(ddiff(i, obs), M1(obs), P1(obs), T1(obs), r1(obs));
                end
                
                %% Setting non-physical values to previous panel values
                set_prev = ~isfinite(M2(i,:)) | ~isreal(M2(i,:));
                
                if any(set_prev)
                    
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
                % latter couldoverwrite flow properties already calculated
                % from other impact methods
                pm = (coni & ~shock) | (obj.shadow(i,:) & any(con(i+1:end, :)));
                
                if any(pm)
                    
                    [Cp2(i, pm), M2(i, pm), P2(i, pm), T2(i, pm), r2(i, pm)] = ...
                        prandtlmeyer(ddiff(i, pm), M1(pm), P1(pm), T1(pm), r1(pm));
                end
            end
            
            % Only return panel results asked for
            obj.Cp = Cp2;
            obj.Mach = M2;
            obj.Pre = P2;
            obj.T = T2;
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
                
                % Anderson p45
                vMp1 = obj.pm_fun(M1);
                % Absolute value as del is relative to flow direction, not
                % convex/concavity
                vMp2 = vMp1 + abs(dTheta);
                
                [M2,~,~] = bisection(obj.pm_fun, M1, 2*maxMach, vMp2);
                                
                [T2_T1, P2_P1, r2_r1] = obj.isentropic(M1, M2);
                
                T2 = T2_T1 .* T1;
                P2 = P2_P1 .* P1;
                r2 = r2_r1 .* r1;
                
                Cp2 = 2 * ((P2/Pinf) - 1)/(g * (Minf^2));
                
                base = isnan(M2);
                
                if any(base)
                    %% TODO: change base function to include this
                    % 0 is default tho?
                    M2(base) = 0;
                    T2(base) = 0;
                    P2(base) = 0;
                    r2(base) = 0;
                    Cp2(base) = obj.base_pressure(Minf, g);
                end
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
            
            obj.Mach(con) = 0;
            obj.T(con) = 0;
            obj.Pre(con) = 0;
            obj.rho(con) = 0;
            obj.Cp(con) = obj.base_pressure(obj.flow.Minf, obj.flow.gamma);
        end
        function [obj, tang, surf] = surface_velocity(obj, unorm, normalise)
            
            tang = crossmat(unorm, obj.flow.U);
            surf = crossmat(tang, unorm);
            
            obj.V = surf;
            
            if nargin >= 3 && normalise
                
                tang_mag = magmat(tang);
                surf_mag = magmat(surf);
                tang = tang./tang_mag;
                surf = surf./surf_mag;
            end
        end
        function Cf = eckert(obj, centre, area, conical)
            
            f = obj.flow;
            g = f.gamma;
            r = f.r;
            
            x = [zeros(1, size(centre, 2)); cumsum(magmat(diff(centre)))];
            
            Te = obj.T;
            Me = obj.Mach;
            % Pe = obj.Pre;
            % Ve = Me * f.a;
            Ve = magmat(obj.V);
            
            %% Eckert
            
            % Assume adiabatic wall T
            % Taw = f.adiabatic(Te);
            % haw = f.cp * Taw;
            h0 = f.cp * f.T0;
            he = f.cp * Te;
            hr = f.cp * f.Tinf + (f.Uinf^2 - Ve.^2)/2 + f.r * (Ve.^2)/2;
            
            Tw = f.Tinf;
            hw = f.cp * Tw;
            tauw = 0.01;
            rw = 8050;
            
            dt = 0.1;
            % Stefan-Boltzman constant
            STF = 5.67e-8;
            % Emissivity
            eps = 0.8;
            F = 1;
            
            for i = 1:10
                
                qw = heating(hw);
                h = qw./(hw - he);
                Tw = Tw + (h .* (hr - hw) - STF * eps * F * Tw.^4)/(rw * f.cp * tauw) * dt;
                hw = Tw * f.cp;
            end
            % Tw = ...
            
            % Original Eckert's reference temperature method
            % Tref = Te * (1 + 0.032*(Me^2) + 0.58*((Tw/Te) - 1));
            
            % Thermoelastic Formulation of a Hypersonic Vehicle Control Surface for Control-Oriented Simulation
            % Tref = Te + 0.5*(Tw - Te) + 0.22*(Tr - Te);
            
            % Meador-Smart reference temperature method - Laminar
            % Tref = Te * (0.45 + 0.55*(Tw/Te) + 0.16*r*((g - 1)/2) * Me^2);
            
            % Meador-Smart reference temperature method - Turbulent
            Tstar = Te .* (0.5*(1 + (Tw./Te)) + 0.16*r*((g - 1)/2) .* Me.^2);
            
            % Ideal gas law
            % Pinf in Anderson2006 for flat plate, Pe?
            rhoStar = f.Pinf./(f.R * Tstar);
            
            % Sutherland's law
            S = 110;
            muStar = f.mu * ((Tstar/f.Tinf).^1.5) .* (f.Tinf + S)./(Tstar + S);
            
            ReStarx = (rhoStar .* Ve .* x)./muStar;
            
            % Accurate up to 10^9 according to above frictionpaper ref
            % cf = (0.37./(log(ReRefx).^2.584)) .* area/Aref;
            
            % Equation from above ref
            % cf = (0.0592./(ReRefx.^0.2)) .* area/Aref;
            
            % Hypersonic and high-temperature gas dynamics (Meador-Smart)
            cf = 0.02296./((ReStarx).^0.139) .* area/obj.Aref;
            
            cf(~isfinite(cf)) = 0;
            

%             StRef = 0.5 * cf * Pr^(-2/3);
%             h = StRef * f.cp * rhoStar * Ve;
            
            % Mangler factor for conical bodies
            if conical
                
                cf = 3^0.5 * cf;
            end
            
            Cf = sum(cf(:));
            
            function qw = heating(hw, stag, laminar)
                
                if nargin < 1 || isempty(hw), hw = 0; end
                if nargin < 2 || isempty(stag), stag = false; end
                if nargin < 3 || isempty(laminar), laminar = true; end
                
                rinf = f.rinf;
                U = f.Uinf;
                d = obj.del;
                
                if stag
                    
                    M = 3;
                    N = 0.5;
                    C = 1.83e-8 * R^-0.5 * (1 - hw/h0);
                    
                elseif laminar
                    
                    M = 3.2;
                    N = 0.5;
                    C = 2.53e-9 * cos(d).^0.5 .* sin(d) .* ...
                        x.^-0.5 .* (1 - hw/h0);
                    
                elseif U > 3962
                    
                    M = 3.7;
                    N = 0.8;
                    C = 2.2e-9 * cos(d).^2.08 .* sin(d).^1.6 .* ...
                        xt.^-0.2 .* (1 - (1.11*hw/h0));
                else
                    M = 3.37;
                    N = 0.8;
                    C = 3.89e-8 * cos(d).^1.78 .* sin(d).^1.6 .* ...
                        xt.^-0.2 .* (Tw/566).^-0.25 * (1 - (1.11*hw/h0));
                end
                
                qw = rinf^N * U^M * C;
            end
        end
%         function stag = stagnation(obj, part)
%             
%             pmax = max(obj.del(1,:));
%             
%             if part.conical
%                 
%                 test = mean(obj.del(1,:));
%             else
%                 row = ceil(size(obj.del, 1));
%                 [stag, id] = max(obj.del(1:row, :));
%             end
%         end

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
            
            for i = 1:numel(fields)
                
                fname = fields{i};
                res = obj.(fname);
                
                if isstring(res) && any(contains(name, res, 'IgnoreCase', true))
                    return
                end
            end
            
            fname = [];
            res = [];
            i = 0;
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
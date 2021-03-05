classdef HASASimp
    
    properties
        
        Ahfp    = 0; % Ratio of horizontal stabiliser area/wing area
        Avfp    = 0; % Ratio of vertical stabiliser area/wing area
        rho_pay = 3.3;
        
        % Wing
        span
        lambda_mid % Mid-chord sweep angle
        t_c % Wing thickness to chord ratio
        lambda % Wing taper ratio
        
        % Propulsion
        Aratio  = 0; % Rocket expansion ratio
        Htsjm   = 0; % Height of scramjet engine
        Nengrt  = 0; % Number of rockets
        Nengtr  = 0; % Number of turboramjets
        Nengrj  = 0; % Number of ramjets
        Nengsj  = 0; % Number of scramjets
        Nengtj  = 0; % Number of turbojets
        Ttotrk  = 0; % Total momentum thrust of all rocket engines (lb)
        Ttott   = 0; % Total momentum thrust of all airbreathing engines (lb)
        Wa      = 0; % Engine airflow (lb/sec)
        
        % Propellant & tank density
        rho_a   = 7; % Vehicle density (Wgtot - Wfuel - Wpay)/Vtot (lb/ft^3)
        rho_f   = 5.25; % Density of hydrogen fuel (lb/ft^3)
        rho_hy  = 0; % Density of hydrazine (lb/ft^3)
        rho_ni  = 0; % Density of nitrogen tetroxide (lb/ft^3)
        rho_rp  = 0; % Density of rp1 (lb/ft^3)
        rho_o   = 0;     % Density of 02 (lb/ft^3)
        rho_th  = 1; % Density of hydrogen tank (lb/ft^3)
        rho_t0  = 0; % Density of oxygen tank (lb/ft^3)
        
        % Misc
        Qmax    = 1000; % Maximum dynamic pressure (lb/ft^2)
        Wins    = 1.5; % Unit weight of thermal protection system (lb/ft^2)
        ULF     = 5; % Ultimate load factor
        mf      = 1.12; % Modifying factor
        del     = 1; % Fuel stored in: fuselage = 1, wings = 0
        vert    = 0; % Vertical take-off
        integral= false; % Integrated fuel tanks
        
        %
        body
        wing
        Sref = 0
        Lb
        Sbtot
        Dbe
        Vtot
        Vpay
        Wpay
        pay_frac
        fuel_frac
        Vfuel
        fuel_ratios = [0 0 0.4 0 0]'
        % Starting point
        mass = 1000;
        
        dynamic_ff = true
        converged = false
    end
    
    properties (Dependent)
        
        AR % Wing aspect ratio
        mass_kg
    end
    
    properties (Constant)
        
        % Usable volume
        nvol = 0.8
        
        m_ft_conv = 3.28084;
        kg_lb_conv = 2.20462;
        lb_kg_conv = 0.45359;
        vol_conv = 0.028317;
        
        ff_table = csvread('HASA_fuelfrac.csv')
    end
    
    methods
        
        function obj = HASASimp(body, wing)
            
            if nargin > 0
                
                obj.body = body;
                
                if nargin > 1
                    
                    obj.wing = wing;
                end
            end
        end
        
        function obj = get_data(obj, b, w)
            
            if nargin < 2 || isempty(b), b = obj.body; end
            if nargin < 3 || isempty(w), w = obj.wing; end
            
            if ~isempty(w)
                
                obj = obj.wing_data(w);
            end

            obj.Dbe = b.width * obj.m_ft_conv;
            obj.Lb = b.length * obj.m_ft_conv;
            obj.Sbtot = 2 * sum(b.quad_data.area(:)) * obj.m_ft_conv^2;
            obj.Vtot = b.volume / obj.vol_conv;
        end
        function obj = wing_data(obj, w)
            
            if nargin < 2 || isempty(w), w = obj.wing; end
            
            % Doubling for half-body
            obj.span = 2 * sum(w.span) * obj.m_ft_conv; % Wing span (ft)
            obj.Sref = 2 * sum(w.area) * obj.m_ft_conv^2;
            
            sections = [w.sections];
            % Already non-dimensionalised (t/c)
            max_t = max([sections.zu] - [sections.zl], [], 1);
            
            % Mean? Weighted mean? See shuttle
            obj.lambda_mid  = mean(w.get_sweep(0.5)); % Mid-chord sweep angle
            obj.t_c         = mean(max_t); % Wing thickness to chord ratio
            obj.lambda      = mean(w.taper); % Wing taper ratio
        end
        function obj = get_payload(obj)
            
            % Useable volume minus fuel volume = volume of payload
            useVol = (obj.nvol * obj.Vtot) - sum(obj.Vfuel);
            % Potential volume, hardcoded maximum (50% of total mass)
            obj.Wpay = min(useVol * obj.rho_pay, 0.5 * obj.mass);
            obj.Vpay = obj.Wpay/obj.rho_pay; % Weight of payload
        end
        function obj = main(obj)
            
            obj.integral = obj.fuel_ratios(4);
            i = 1;
            
            while true
                
                old = obj.mass;
                [obj, new] = obj.weigh();
                
                % Convergence/penalty criteria
                conv = [abs(new - old) < 1, i > 1000, ~isfinite(new)];
                obj.mass = new;
                
                if any(conv)
                    
                    obj.converged = conv(1);
                    
                    if ~obj.converged 
                        
                        % error('weight estimation did not converge')
                        % Apply penalty if i over number of allowed iters
                        obj.mass = inf;
                        obj.Vfuel = inf;
                    else
                        % Max volume = usable percentage * total
                        maxVol = HASASimp.nvol * obj.Vtot;

                        payVol = obj.Vpay;
                        fuelVol = sum(obj.Vfuel);
                        
                        if payVol + fuelVol > maxVol
                            
                            % Reset payload volume if needed
                            payVol = maxVol - fuelVol;
                            obj.Vpay = payVol;
                            obj.Wpay = obj.Vpay * obj.rho_pay;
                        end
                    end
                    
                    break
                else
                    i = i + 1;
                end
            end
            
            obj.pay_frac = obj.Wpay/obj.mass;
        end
        function [obj, Wgtot] = weigh(obj)
            
            Wgtot = obj.mass;
            
            %% TODO: Need values
            rho_t = [0, 0, obj.rho_th, obj.rho_t0, 0]'; % ?
            rho_fuel = [obj.rho_hy, obj.rho_ni, obj.rho_f, obj.rho_o, obj.rho_rp]';
            
            %%
            if obj.dynamic_ff
                %% TODO: Getting negative values for low weight, HASA graph look positive even at 0
                % Minimum of 0.1 fuel fraction
                ff_tot = interp1(obj.ff_table(:,2), ...
                    obj.ff_table(:,1), Wgtot, 'linear', 'extrap');
                
                % Fuel fraction limited within reasonable bounds
                ff_tot = min(max(ff_tot, 0.2), 0.8);
                
                % Normalise fuel ratios and multiplying by total fuel frac
                % Thus sum(ff) == ff_tot
                ff = (obj.fuel_ratios/sum(obj.fuel_ratios)) * ff_tot;
                
                obj.fuel_frac = ff_tot;
            else
                ff = obj.fuel_ratios;
                obj.fuel_frac = sum(ff);
            end
            
            Wfuel = Wgtot * ff;
            
            obj.Vfuel = Wfuel./rho_fuel;
            % Tank density and volume of fuel
            Wtnk = sum(rho_t .* obj.Vfuel);
            
            %% Payload defined by fraction of total mass, or by body volume
            if ~isempty(obj.pay_frac)
                
                obj.Wpay = Wgtot * obj.pay_frac;
                obj.Vpay = obj.Wpay/obj.rho_pay;
            else
                obj = obj.get_payload;
            end
            
            %% Structural
            
            Wemp = Wgtot - sum(Wfuel);
            
            if obj.integral
                
                Wemp = max(Wemp - Wtnk, 1);
                Wb = Wtnk;
                Wtnk = 0;
            else
                o = ((obj.Lb * obj.ULF)/obj.Dbe)^0.15 * obj.Qmax^0.16 * ...
                    obj.Sbtot^1.05;
                Wb = 0.341 * obj.mf * o;
            end
            
            Swfh = obj.Sref * obj.Ahfp; % Horizontal stabilizer planform area (ft^2)
            Swfv = obj.Sref * obj.Avfp; % Vertical stabiliser planform area (ft^2)
            
            if obj.Sref
                
                Ww = 0.2958 * obj.mf * (((Wemp * obj.ULF)/1000)^0.52 ...
                    * obj.Sref^0.7 * obj.AR^0.47 * ((1 + obj.lambda)/(obj.t_c))^0.4 ...
                    * (0.3 + (0.7/cos(obj.lambda_mid))))^1.017;
          
                Lam = (Wgtot/obj.Sref)^0.6 * Swfh^1.2 * obj.Qmax^0.8;
                Wth = 0.00035 * Lam;
                Wtv = 5 * Swfv^1.09;
            else
                %% TODO: 
                % Future release will need way of determining tail weights
                % on lifting bodies
                Ww = 0;
                Wth = 0;
                Wtv = 0;
            end
            
            Wtps = obj.Wins * (obj.Sbtot/2 + obj.Sref + Swfh);
            
            if obj.vert
                
                Wg = Wemp;
            else
                Wg = Wgtot;
            end
            
            Wgear = 0.00916 * Wg^1.124;
            
            % [Airbreathing, Rocket]
            Wthrst = [0.00625 * obj.Ttott + 69, 0.0025 * obj.Ttotrk];
            
            Wstr = Wb + Ww + Wth + Wtv + Wtps + Wgear + sum(Wthrst);
            
            %% Engines
            Wttj = obj.Nengtj * (obj.Wa * 133.3 - 16600)/4;
            Wttr = obj.Nengtr * 1782.63 * exp(0.003*obj.Wa);
            Wtrj = 0.01 * obj.Ttott; % * Nengrj;
            Wtsj = obj.Nengsj * (87.5 * obj.Htsjm - 850);
            Wtrt = 0.00766 * obj.Ttotrk + 0.00033 * obj.Ttotrk ...
                * obj.Aratio^0.5 + 130 * obj.Nengrt;
            
            Wpros = Wttj + Wttr + Wtrj + Wtsj + Wtrt + Wtnk;
            
            %% Hydraulics/Avionics/Electronics
            if obj.Sref
                
                psi = (((obj.Sref + Swfv + Swfh) * obj.Qmax)/1000)^0.334 ...
                    * (obj.Lb + obj.span)^0.5;
            else
                psi = (((obj.Sbtot + Swfv + Swfh) * obj.Qmax)/1000)^0.334 ...
                    * obj.Lb^0.5;
            end
            
            Whydr = 2.64 * psi;
            Wtavcs = 66.37 * Wgtot^0.361;
            
            O = Wgtot^0.5 * obj.Lb^0.25;
            Welect = 1.167*O;
            
            Wequip = 1e4 + 0.01*(Wgtot - 3e-7);
            
            Wsub = Whydr + Wtavcs + Welect + Wequip;
            
            Wgtot = obj.Wpay + sum(Wfuel) + Wstr + Wpros + Wsub;
        end
        function a = check(obj)
            %% Feasibility checks
            % fuel fraction not too high, payload mass > 0
            a = [obj.fuel_frac, obj.Wpay];
        end
        function obj = set.Vfuel(obj, in)
            
            % Can be nan due to div by 0 for unknown fuel densities
            in(isnan(in)) = 0;
            obj.Vfuel = in;
        end
        function a = get.AR(obj)
            
            if obj.Sref
                
                a = obj.span^2/obj.Sref;
            else
                a = 0;
            end
        end
        function a = get.mass_kg(obj)
            
            a = obj.mass * obj.lb_kg_conv;
        end
    end
    
    methods (Static)
        
        function [init_obj, a] = define()
            
            a(1) = struct('name', "pay_frac", 'min', 0.02, 'max', 0.5);
            
            init_obj = HASASimp();
            [init_obj, a] = Geometry.define(init_obj, a);
        end
        function validate()
            
            table1 = xlsread('HASAdata', 1);
            table2 = xlsread('HASAdata', 2);
            
            for i = 1:size(table2, 2)
                
                data1 = table1(:,i);
                data2 = table2(:,i);
                
                obj = HASASimp();
                
                % Other inputs? (Table 2)
                obj.fuel_ratios   = data1(1:5);
                obj.span        = data2(3); % Wing span (ft)
                obj.Sref        = data2(5); % Wing reference area (ft^2)
                
                % Input list (Table 3)
                % Geometry
                obj.Ahfp        = data1(6); % Ratio of horizontal stabiliser area/wing area
                obj.Avfp        = data1(7); % Ratio of vertical stabiliser area/wing area
                obj.lambda_mid  = deg2rad(data1(12)); % Mid-chord sweep angle
                obj.t_c         = data1(15); % Wing thickness to chord ratio
                obj.lambda      = data1(16); % Wing taper ratio
                obj.Vpay        = data1(17); % Volume of payload
                obj.Wpay        = data1(18); % Weight of payload
                
                % Propulsion
                obj.Aratio      = data1(20); % Rocket expansion ratio
                obj.Htsjm       = data1(21); % Height of scramjet engine
                obj.Nengrt      = data1(22); % Number of rockets
                obj.Nengtr      = 0; % Number of turboramjets
                obj.Nengrj      = 0; % Number of ramjets
                obj.Nengsj      = data1(23); % Number of scramjets
                obj.Nengtj      = data1(24); % Number of turbojets
                obj.Ttotrk      = data1(25); % Total momentum thrust of all rocket engines (lb)
                obj.Ttott       = data1(26); % Total momentum thrust of all airbreathing engines (lb)
                obj.Wa          = data1(27); % Engine airflow (lb/sec)
                
                % Propellant & tank density
                obj.rho_a       = data1(28); % Vehicle density (Wgtot - Wfuel - Wpay)/Vtot (lb/ft^3)
                obj.rho_f       = data1(29); % Density of hydrogen fuel (lb/ft^3)
                obj.rho_hy      = data1(30); % Density of hydrazine (lb/ft^3)
                obj.rho_ni      = data1(31); % Density of nitrogen tetroxide (lb/ft^3)
                obj.rho_rp      = data1(32); % Density of rp1 (lb/ft^3)
                obj.rho_o       = 71.27;     % Density of 02 (lb/ft^3)
                obj.rho_th      = data1(33); % Density of hydrogen tank (lb/ft^3)
                obj.rho_t0      = data1(34); % Density of oxygen tank (lb/ft^3)
                
                % Misc
                obj.Qmax        = data1(35); % Maximum dynamic pressure (lb/ft^2)
                obj.Wins        = data1(36); % Unit weight of thermal protection system (lb/ft^2)
                obj.ULF         = data1(37); % Ultimate load factor
                obj.mf          = data1(38); % Modifying factor
                obj.del         = data1(39); % Fuel stored in: fuselage = 1, wings = 0
                obj.vert        = data1(40); % Vertical take-off
                
                % Additional/unknown inputs
                obj.mass        = 1000; % INPUT
                obj.Sbtot       = data2(13); % INPUT
                obj.Lb          = data2(2);
                obj.Dbe         = data2(4);
                
                obj = obj.main();
                
                save(i) = obj;
            end
        end
        function obj = test()
            
            obj = HASASimp();
            obj.Sref = 357.5; % ft2
            obj.span = 27.7083; % ft
            obj.lambda = 0.11;
            obj.lambda_mid = deg2rad(45);
            obj.t_c = 0.1;
            obj.Sbtot = 8000;
            
%             vol = 4*obj.Sbtot;
%             obj = obj.get_payload(vol);
            
            obj.pay_frac = 0.05;
            obj.Lb = 17.57 * HASASimp.m_ft_conv;
            obj.Dbe = 1.82 * HASASimp.m_ft_conv;
            
            obj = obj.main();
        end
    end
end
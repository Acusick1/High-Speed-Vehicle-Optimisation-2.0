classdef Parsec < Aerofoil
    
    properties
        
        rleu
        xup
        zup
        zxxup
        xlo
        zlo
        zxxlo
        zte
        dzte
        ate
        bte
        rlel
    end
    
    methods
        
        function obj = Parsec(p)
            
            if nargin >= 1
                % Coefficient Matrix Upper
                Cup = [1 1 1 1 1 1;
                    p(2)^(1/2) p(2)^(3/2) p(2)^(5/2) p(2)^(7/2) p(2)^(9/2) p(2)^(11/2);
                    1/2 3/2 5/2 7/2 9/2 11/2;
                    1/2*p(2)^(-1/2) 3/2*p(2)^(1/2) 5/2*p(2)^(3/2) 7/2*p(2)^(5/2) 9/2*p(2)^(7/2) 11/2*p(2)^(9/2);
                    -1/4*p(2)^(-3/2) 3/4*p(2)^(-1/2) 15/4*p(2)^(1/2) 35/4*p(2)^(3/2) 63/4*p(2)^(5/2) 99/4*p(2)^(7/2);
                    1 0 0 0 0 0];
                
                % Coefficient Matrix Lower
                Clo = [1 1 1 1 1 1;
                    p(5)^(1/2) p(5)^(3/2) p(5)^(5/2) p(5)^(7/2) p(5)^(9/2) p(5)^(11/2);
                    1/2 3/2 5/2 7/2 9/2 11/2;
                    1/2*p(5)^(-1/2) 3/2*p(5)^(1/2) 5/2*p(5)^(3/2) 7/2*p(5)^(5/2) 9/2*p(5)^(7/2) 11/2*p(5)^(9/2);
                    -1/4*p(5)^(-3/2) 3/4*p(5)^(-1/2) 15/4*p(5)^(1/2) 35/4*p(5)^(3/2) 63/4*p(5)^(5/2) 99/4*p(5)^(7/2);
                    1 0 0 0 0 0];
                
                % Coefficients of 'b' upper and lower
                Bup = [p(8)+p(9)/2 p(3) tand((2*p(10)-p(11))/2) 0 p(4) sqrt(2*p(1))]';
                
                Blo = [p(8)-p(9)/2 p(6) tand((2*p(10)+p(11))/2) 0 p(7) -sqrt(2*p(12))]';
                
                Aup = Cup\Bup;
                Alo = Clo\Blo;
                
                % Parsec AnAloyticAlo Formulation
                zu = zeros(size(obj.xu));
                zl = zeros(size(obj.xl));
                
                for k = 1:6
                    
                    zu = zu + Aup(k) .* obj.xu.^(k - 0.5);
                    zl = zl + Alo(k) .* obj.xl.^(k - 0.5);
                end
                
                obj.zu = zu;
                obj.zl = zl;
            end
        end
        function obj = generate(obj)
            
            p(1) = obj.rleu;
            p(2) = obj.xup;
            p(3) = obj.zup;
            p(4) = obj.zxxup;
            p(5) = obj.xlo;
            p(6) = obj.zlo;
            p(7) = obj.zxxlo;
            p(8) = obj.zte;
            p(9) = obj.dzte;
            p(10) = obj.ate;
            p(11) = obj.bte;
            p(12) = obj.rlel;
            
            % Coefficient Matrix Upper
            Cup = [1 1 1 1 1 1;
                p(2)^(1/2) p(2)^(3/2) p(2)^(5/2) p(2)^(7/2) p(2)^(9/2) p(2)^(11/2);
                1/2 3/2 5/2 7/2 9/2 11/2;
                1/2*p(2)^(-1/2) 3/2*p(2)^(1/2) 5/2*p(2)^(3/2) 7/2*p(2)^(5/2) 9/2*p(2)^(7/2) 11/2*p(2)^(9/2);
                -1/4*p(2)^(-3/2) 3/4*p(2)^(-1/2) 15/4*p(2)^(1/2) 35/4*p(2)^(3/2) 63/4*p(2)^(5/2) 99/4*p(2)^(7/2);
                1 0 0 0 0 0];
            
            % Coefficient Matrix Lower
            Clo = [1 1 1 1 1 1;
                p(5)^(1/2) p(5)^(3/2) p(5)^(5/2) p(5)^(7/2) p(5)^(9/2) p(5)^(11/2);
                1/2 3/2 5/2 7/2 9/2 11/2;
                1/2*p(5)^(-1/2) 3/2*p(5)^(1/2) 5/2*p(5)^(3/2) 7/2*p(5)^(5/2) 9/2*p(5)^(7/2) 11/2*p(5)^(9/2);
                -1/4*p(5)^(-3/2) 3/4*p(5)^(-1/2) 15/4*p(5)^(1/2) 35/4*p(5)^(3/2) 63/4*p(5)^(5/2) 99/4*p(5)^(7/2);
                1 0 0 0 0 0];
            
            %Coefficients of 'b' upper and lower
            Bup = [p(8)+p(9)/2 p(3) tand((2*p(10)-p(11))/2) 0 p(4) sqrt(2*p(1))]';
            
            Blo = [p(8)-p(9)/2 p(6) tand((2*p(10)+p(11))/2) 0 p(7) -sqrt(2*p(12))]';
            
            Aup = Cup\Bup;
            Alo = Clo\Blo;
            
            % Parsec AnAloyticAlo Formulation
            zu = zeros(size(obj.xu));
            zl = zeros(size(obj.xl));
            
            for k = 1:6
                
                zu = zu + Aup(k) .* obj.xu.^(k - 0.5);
                zl = zl + Alo(k) .* obj.xl.^(k - 0.5);
            end
            
            obj.zu = zu;
            obj.zl = zl;
        end
        function obj = set.rlel(obj, val)
            
            obj.rlel = val;
            obj = obj.generate;
        end
    end
    
    methods (Static)
        
        function test()
            %% NACA2411 Test Parameters
            comp = Aerofoil.getaerofoil("NACA2411.dat");
            comp = Aerofoil(comp);
            
            rleu = 0.0216;              % Rleu: Upper Leadng edge Radius(upper)
            xup = 0.3445;               % Xup: upper crest position in horizontAlo coordinates
            zup = 0.07912;              % Zup: upper crest position in verticAlo coordinates
            zxxup = -0.6448;              % ZXXup: Upper Crest Curvature
            xlo = 0.17;               % Xlo: lower crest position in horizontAlo coordinates
            zlo = -0.033797;              % Zlo: lower crest position in verticAlo coordinates
            zxxlo = 0.6748;               % ZXXlo: Lower Crest Curvature
            zte = 0;
            dzte = 0;
            ate = -4.785;             % P10: Trailing Edge direction angle
            bte = 15.082;                % P11: TrAlong Edge Wedge Angle
            rlel = 0.008;             % Rlel: Lower Leading edge radius (lower)
            
            %% Application of swarm approach and artificial neural networks for airfoil shape optimization
            % Unsure of specific aerofoils
            % rleu = 0.01;
            % xup = 0.3;
            % zup = 0.06;
            % zxxup = -0.45;
            % xlo = 0.3;
            % zlo = -0.06;
            % zxxlo = 0.45;
            % zte = 0;
            % dzte = 0;
            % ate = 0;
            % bte = 17;
            % rlel = 0.01;
            
            % rleu = 0.01;
            % xup = 0.4324;
            % zup = 0.063;
            % zxxup = -0.4363;
            % xlo = 0.3438;
            % zlo = -0.0589;
            % zxxlo = 0.7;
            % zte = 0;
            % dzte = 0;
            % ate = -6.8;
            % bte = 8.07;
            % rlel = 0.01;
            
            %%
            p(1) = rleu;
            p(2) = xup;
            p(3) = zup;
            p(4) = zxxup;
            p(5) = xlo;
            p(6) = zlo;
            p(7) = zxxlo;
            p(8) = zte;
            p(9) = dzte;
            p(10) = ate;
            p(11) = bte;
            p(12) = rlel;
            
            obj = Parsec(p);
            a = obj.checks;
            obj.plotter(obj, comp);
            
        end
        function [init, a] = define()
            % Standard Parsec optimisation definition
            
            % var_min = [0.001, 0.2, 0.02, -1.2, 0.2, -0.08, 0, -0.02, 0, -25, 3, 0.001];
            % var_max = [0.1, 0.7, 0.12, 0, 0.7, -0.02, 1.2, 0.02, 0.02, 2, 40, 0.1];
            var_min = [0.001, 0.05, 0.005, -1.2, 0.05, -0.1, 0, -0.02, 0, -25, 3, 0.001];
            var_max = [0.1, 0.8, 0.2, 0, 0.8, -0.02, 1.2, 0.02, 0.02, 2, 40, 0.1];
            name = ["rleu", "xup", "zup", "zxxup", "xlo", "zlo",...
                "zxxlo", "zte", "dzte", "ate", "bte", "rlel"];
            con = [];
            trans = [];
            opt = true(size(name));
            
            init = Parsec();
            a.name = name(:)';
            a.min = var_min(:)';
            a.max = var_max(:)';
            a.val = (var_min(:) + var_max(:))'/2;
            % a = OptVariables(var_min, var_max, name, con, trans, opt);
        end
    end
end
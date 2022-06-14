classdef Forces < Combinable
    
    properties
    
        L   % Lift (N)
        D   % Drag (N)
        Df  % Friction drag (N)
        M   % Pitching moment (Nm) 
    end
    
    methods
        function self = Forces(force_coeffs, q, Aref)
            %FORCES constructor
            %   Inputs:
            %   force_coeffs - single or array of ForceCoeffs objects, with
            %       arrays defining multiple objects part of a single
            %       configuration that will be combined to produce overall
            %       forces
            %   q - Dynamic pressure (0.5 * rho * U^2)
            %   Aref - reference area
            
            if nargin > 0
                
                self.L  = q * sum([force_coeffs.Cl]) * Aref;
                self.D  = q * sum([force_coeffs.Cd]) * Aref;
                self.Df = q * sum([force_coeffs.Cf]) * Aref;
                self.M  = q * sum([force_coeffs.Cm]) * Aref;
            end
        end
    end
end
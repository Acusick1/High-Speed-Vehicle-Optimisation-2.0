classdef ForceCoeffs < Combinable
    
    properties
        
        Cl      % Lift coefficient
        Cdp     % Inviscid pressure drag coefficient
        Cm      % Pitching moment coefficient
        CN      % Normal coefficient
        CA      % Axial coefficient
        Cf      % Friction drag coefficient
        Cd      % Overall drag coefficient
    end
    
    methods
        
        function obj = ForceCoeffs(Cp, Cf, area, unit_norm, alpha, Aref)
            %FORCECOEFFS constructor
            %   Inputs:
            %   Cp - inviscid pressure coefficient per panel [NxM]
            %   Cf - viscous pressure coefficient per panel [NxM]
            %   area - panel areas [NxM]
            %   unit_norm - panel unit normals [NxMx3]
            %   alpha - angle of attack (rad)
            %   Aref - reference area (m^2)
            if nargin > 0
                
                % Angles between body and inertial NED axes
                xyAngle = alpha;
                xzAngle = 0;
                yzAngle = pi/2 - alpha;

                nx = unit_norm(:,:,1);
                ny = unit_norm(:,:,2);
                nz = unit_norm(:,:,3);
                
                %% Create temporary unnormalised object
                % Aerodynamic characteristics
                % Development of an Aerodynamics Code for the Optimisation of Hypersonic Vehicles
                f.Cl = -((Cp .* area .* nx) * sin(xyAngle)) -...
                    ((Cp .* area .* ny) * sin(xzAngle)) -...
                    ((Cp .* area .* nz) * sin(yzAngle));
                
                f.Cdp = -((Cp .* area .* nx) * sin(yzAngle)) +...
                    ((Cp .* area .* ny) * sin(xzAngle)) -...
                    ((Cp .* area .* nz) * sin(xyAngle));
                
                if isempty(Cf)
                    f.Cd = f.Cdp;
                    f.CA = -Cp .* area .* nx;
                else
                    f.Cf = Cf .* area;
                    f.Cd = f.Cdp + f.Cf;
                    f.CA = -(Cp + Cf) .* area .* nx;
                end
                
                f.CN = -Cp .* area .* nz;
                f.Cm = f.CA - f.CN;
                
                %% Sum and normalise coefficients to ForceCoeffs object
                fn = fieldnames(f);
                
                for i = 1:numel(fn)
                    
                    obj.(fn{i}) = sum(f.(fn{i})(:))/Aref;
                end
            end
        end
    end
end
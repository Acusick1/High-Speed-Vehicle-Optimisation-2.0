classdef ForceCoeffs < Combinable
    
    properties
        
        Cl
        Cdp
        Cm
        CN
        CA
        Cf
        Cd
        Cl_Cd
    end
    
    methods
        
        function obj = ForceCoeffs(data, Cp, Cf, alpha, Aref)
            
            if nargin > 0
                
                xyAngle = alpha;
                xzAngle = 0;
                yzAngle = pi/2 - alpha;
                
                area = data.area;
                unit_norm = data.unit_norm;
                
                %%
                nx = unit_norm(:,:,1);
                ny = unit_norm(:,:,2);
                nz = unit_norm(:,:,3);
                
                % Part aerodynamic characteristics
                % Development of an Aerodynamics Code for the Optimisation of
                % Hypersonic Vehicles
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
                
                % f.Cm = -(Cp .* area .* nx) + (Cp .* area .* nz);
                f.CN = -Cp .* area .* nz;
                % f.CA = -Cp .* area .* nx;
                f.Cm = f.CA - f.CN;
                
                %% Sum and normalise coefficients
                fn = fieldnames(f);
                
                for i = 1:numel(fn)
                    
                    obj.(fn{i}) = sum(f.(fn{i})(:))/Aref;
                end
            end
        end
    end
end
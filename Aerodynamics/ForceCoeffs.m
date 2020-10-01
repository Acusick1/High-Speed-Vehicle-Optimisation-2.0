classdef ForceCoeffs < Combinable
    
    properties
        
        Cl
        Cd
        Cm
        CN
        CA
        Cf
        Aref
    end
    
    methods
        
        function obj = ForceCoeffs(data, Cp, Cf, alpha, Aref)
            
            if nargin > 0
                
                xyAngle = alpha;
                xzAngle = 0;
                yzAngle = pi/2 - alpha;
                
                area = data.area;
                unit_norm = data.unitNorm;
                
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
                
                f.Cd = -((Cp .* area .* nx) * sin(yzAngle)) +...
                    ((Cp .* area .* ny) * sin(xzAngle)) -...
                    ((Cp .* area .* nz) * sin(xyAngle));
                
                f.Cm = -(Cp .* area .* nx) + (Cp .* area .* nz);
                f.CN = -Cp .* area .* nz;
                f.CA = -Cp .* area .* nx;
                f.Cf = Cf;
                
                %% Sum and normalise coefficients
                fn = fieldnames(f);
                
                for i = 1:numel(fn)
                    
                    obj.(fn{i}) = sum(f.(fn{i})(:))/Aref;
                end
                
                obj.Aref = Aref;
            end
        end
    end
end
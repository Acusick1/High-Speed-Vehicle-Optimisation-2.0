classdef Forces < Combinable
    
    properties
    
        L
        D
        M
    end
    
    methods
        function obj = Forces(fco, flow)
            
            if nargin > 0
                
                rinf = flow.rinf;
                Uinf = flow.Uinf;
                
                obj.L = rinf * (Uinf^2) * sum([fco.Cl]) * fco.Aref;
                obj.D = rinf * (Uinf^2) * sum([fco.Cd]) * fco.Aref;
                obj.M = rinf * (Uinf^2) * sum([fco.Cm]) * fco.Aref;
            end
        end
    end
end
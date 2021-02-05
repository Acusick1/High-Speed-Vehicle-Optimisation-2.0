classdef Forces < Combinable
    
    properties
    
        L
        D
        M
    end
    
    methods
        function obj = Forces(fco, flow, Aref)
            
            if nargin > 0
                
                rinf = flow.rinf;
                Uinf = flow.Uinf;
                
                obj.L = rinf * (Uinf^2) * sum([fco.Cl]) * Aref;
                obj.D = rinf * (Uinf^2) * sum([fco.Cd]) * Aref;
                obj.M = rinf * (Uinf^2) * sum([fco.Cm]) * Aref;
            end
        end
    end
end
classdef Forces < Combinable
    
    properties
    
        L
        D
        Df
        M
    end
    
    methods
        function obj = Forces(fco, flow, Aref)
            
            if nargin > 0
                
                rinf = flow.rinf;
                Uinf = flow.Uinf;
                
                obj.L = 0.5 * rinf * (Uinf^2) * sum([fco.Cl]) * Aref;
                obj.D = 0.5 * rinf * (Uinf^2) * sum([fco.Cd]) * Aref;
                obj.Df = 0.5 * rinf * (Uinf^2) * sum([fco.Cf]) * Aref;
                obj.M = 0.5 * rinf * (Uinf^2) * sum([fco.Cm]) * Aref;
            end
        end
    end
end
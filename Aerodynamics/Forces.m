classdef Forces < Combinable
    
    properties
    
        L
        D
        Df
        M
    end
    
    methods
        function self = Forces(fco, flow, Aref)
            
            if nargin > 0
                
                rinf = flow.rinf;
                Uinf = flow.Uinf;
                
                self.L = 0.5 * rinf * (Uinf^2) * sum([fco.Cl]) * Aref;
                self.D = 0.5 * rinf * (Uinf^2) * sum([fco.Cd]) * Aref;
                self.Df = 0.5 * rinf * (Uinf^2) * sum([fco.Cf]) * Aref;
                self.M = 0.5 * rinf * (Uinf^2) * sum([fco.Cm]) * Aref;
            end
        end
    end
end
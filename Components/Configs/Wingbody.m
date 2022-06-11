classdef Wingbody < Configuration
    %WINGBODY is a metaclass configuration this merges wing and body
    %components to a unified wingbody. Optionally, a system component can 
    %be added to compute vehicle characteristics.
    properties
        
        % Outputs
        lower_body
        upper_body
        outer_wing
        full_wing
        Aref
        mass_kg
        pay_frac
    end
    
    methods
        function self = generate(self)
            
            body = self.parts{self.type == "Body"};
            wing = self.parts{self.type == "Wing"};
            sys = self.parts{self.type == "HASA"};
            
            %% TODO: Just pass aerofoil to wing as cell?
            aerofoil_cell = self.parts(self.type == "Aerofoil");
            for i = numel(aerofoil_cell):-1:1
                
                aerofoil(i) = aerofoil_cell{i}.generate();
            end
            
            wing.sections = aerofoil;
            body = body.generate();
            
            [self.outer_wing, self.lower_body, self.upper_body, self.full_wing] = ...
                vec_wingbody(wing, body);
            
            if ~isempty(sys)
                %% TODO Don't attach body and wing, just run program
                % Attach old (full) body and updated wing
                sys.body = body;
                sys.wing = self.outer_wing;
                sys = sys.get_data;
                sys = sys.main;
                
                self.mass_kg = sys.mass_kg;
                self.pay_frac = sys.pay_frac;
            end
            
            if isempty(self.outer_wing)
                
                self.Aref = body.area;
            else
                self.Aref = sum(self.full_wing.area);
            end

        end
    end
end
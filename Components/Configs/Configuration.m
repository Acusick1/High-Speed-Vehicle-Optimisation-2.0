classdef Configuration
    %CONFIGURATION is a metaclass that takes a number of geometric/system
    %based components to define a unified configuration.
    %   Configurations must define their own generate method, since
    %   generating individual components is not enough to create a unified
    %   definition. The generate method contains the necessary logic to
    %   combine individual components.
    properties
        
        % Inputs
        parts
        type
    end
    
    methods
        function self = Configuration(parts, type)
            %CONFIGURATION takes the necessary parts along with their
            %corresponding types
            %   Inputs:
            %   parts - cell of objects that makeup configuration
            %   type - defines each part type so that they can be unpacked
            %       in generation method
            if nargin > 0
                
                self.parts = parts;
                self.type = type;
            end
        end
        
        function generate(self)
            error("Generate method not implemented for subclass: %s", class(self))
        end
    end
end